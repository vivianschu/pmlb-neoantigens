# How Somatic Mutations Are Integrated into the Neoantigen Pipeline

This document walks through, step by step, how a somatic variant call (one row
in a MAF file) becomes a mutant protein context in the FragPipe search database.

The relevant code lives in [workflow/mutation_reference/](.). The user-facing
entry points are [01_build_mutated_db.sh](01_build_mutated_db.sh) (SLURM
wrapper) and [build_mutated_protein_db.py](build_mutated_protein_db.py) (the
Python CLI). Everything described below happens inside Step 01.

---

## 1. The biological problem

Mass-spectrometry immunopeptidomics identifies peptides eluted from MHC. To
detect *neoantigens*, the database used by the search engine must contain the
mutant peptide sequences a tumour can present. A standard reference proteome
(e.g. UniProt) will not contain them — by definition, neoantigens differ from
the germline.

The pipeline therefore builds a **sample-specific, variant-centred mutant
protein FASTA**: for every coding somatic mutation in a tumour, the local
mutant protein context is added to the search space. Only the altered region
plus a flank for peptide-length coverage is written; the rest of the protein
is unchanged and already covered by the reference proteome, so writing it
again would only inflate the search space.

---

## 2. Inputs

| Input | Source | Role |
|---|---|---|
| Somatic MAF | `Org_exome.data_mutations_extended.gt4.202507.txt` | One row per somatic variant call, VEP-annotated, FILTER=PASS or common_variant |
| cDNA ORF FASTA | Output of [00_fetch_cdna_orfs.sh](00_fetch_cdna_orfs.sh) — derived from GENCODE v47 `pc_transcripts` | Reference CDS in nucleotide space, headers `>GENE|TRANSCRIPT_ID|ORF` |
| Sample manifest (optional) | `data/rna/260313_manifest_final.tsv` | Restricts processing to a subset of samples |
| HLA table (optional) | `data/rna/hla_summary.csv` (arcasHLA output) | Drives peptide×HLA pairing in the ML training table |
| TPM / RNA VAF tables (optional) | Bulk RNA-seq outputs | Filter or annotate variants by expression / RNA support |

**Why cDNA, not protein?** HGVSc coordinates (`c.743G>A`) are defined relative
to the start of the CDS. Applying mutations at the cDNA level lets us handle
SNVs, MNVs, in-frame indels, and frameshifts in one uniform framework and
then translate the result — instead of trying to special-case each variant
type at the protein level.

---

## 3. Building the cDNA ORF reference (Step 00, run once)

[00_fetch_cdna_orfs.sh](00_fetch_cdna_orfs.sh) takes the GENCODE v47
protein-coding-transcript FASTA — which contains full mRNAs
(5′ UTR + CDS + 3′ UTR) — and extracts just the CDS for each transcript:

1. Parse the `CDS:start-end` field in the GENCODE header (1-based inclusive).
2. Slice the transcript sequence to those coordinates.
3. If the resulting sequence does not start with `ATG`, fall back to the first
   in-transcript ATG. If neither works, skip the transcript.
4. Write `>GENE|TRANSCRIPT_ID|ORF` headers to
   `mutation_reference/orf/cdna_orfs.fasta`.

The result is a clean, ATG-anchored cDNA ORF reference whose coordinates
align with HGVSc.

[orf.py](orf.py) loads this FASTA at runtime, validates each record
(divisibility by 3, internal stops, invalid bases), and indexes records by
both transcript ID and gene symbol.

---

## 4. Loading and pre-filtering the MAF

[mutations.py](mutations.py) reads the MAF with `pandas`. For each row:

- **Sample manifest gating** — if a manifest was passed, rows whose
  `Tumor_Sample_Barcode` is not in it are dropped. Suffixes such as `_pmlb`
  or `_novo` are stripped before matching so RNA-seq sample IDs collapse
  onto MAF barcodes.
- **Silent classifications dropped** — `Silent`, `Synonymous`,
  `Synonymous_Mutation`, `RNA` are removed before any protein-level work
  begins. They cannot generate neoantigens.
- **HGVSc required** — rows without a parsable HGVSc string are skipped
  (the pipeline operates entirely in CDS coordinates).

Each retained row becomes one `Mutation` dataclass instance carrying the
sample, gene, transcript, HGVSc/HGVSp, variant type, and (where available)
tumour VAF, alt-read count, and depth.

---

## 5. Parsing HGVSc into coordinates

[`parse_hgvsc`](mutations.py) in [mutations.py](mutations.py) is a regex
parser handling the variant grammars that appear in this MAF:

| HGVSc | Type | Parsed as `(start, end, ref, alt, type)` |
|---|---|---|
| `c.743G>A` | SNV | `(743, 743, "G", "A", "SNV")` |
| `c.743_744GC>TA` | MNV | `(743, 744, "GC", "TA", "MNV")` |
| `c.743delG` / `c.743_745delGCC` | DEL | `(743, 745, "GCC", "", "DEL")` |
| `c.743_744insAA` | INS | `(743, 744, "", "AA", "INS")` |
| `c.743delinsAA` | INDEL | `(743, 743, "", "AA", "INDEL")` |
| `c.743dupA` | DUP → INS | `(743, 743, "", "A", "INS")` |

**The `N` placeholder.** This MAF often writes `c.157delNNNNN` or
`c.10138N>A` because the upstream caller did not record the explicit
reference base. The parser detects an all-`N` ref and clears it. Later, the
true reference is read from the ORF sequence at that position — so the call
is still applied correctly, just without a redundant cross-check.

UTR (`c.*`, `c.-`) and intronic positions are filtered here.

---

## 6. Mapping a variant to its ORF

For each `(sample, gene, transcript)` group, [`_find_orf`](build_mutated_protein_db.py)
selects the ORF record to mutate:

1. Look up the MAF `Transcript_ID` in the transcript index (with and without
   version suffix).
2. If not found, fall back to the gene symbol and pick the **longest** ORF
   among that gene's transcripts. This is a deliberate choice: peptide
   evidence is identification-by-sequence, so picking the longest CDS
   maximises the chance that any presented peptide is representable.
3. If neither match succeeds, the variant is logged with
   `status=skipped_no_matching_orf` and dropped.

---

## 7. Optional RNA / expression filtering

If TPM and/or RNA variant tables are loaded, each variant is checked against
the user-supplied thresholds (`--min-tpm`, `--min-rna-vaf`, `--min-rna-depth`).
By default, **missing** RNA data does not fail the variant — the variant
passes through with blank RNA columns. Set `--require-rna-filters` to flip
this and treat missing data as a filter failure.

Filter outcomes are written per-variant to `rna_feature_table.tsv` regardless
of whether they pass, so failures can be inspected post hoc.

---

## 8. Applying the variant to the cDNA

[`apply_mutations_to_sequence`](mutations.py) takes the ORF sequence and the
list of `Mutation` objects belonging to that ORF and produces a mutated
nucleotide string.

Key mechanics:

- The sequence is materialised as a list of characters so insertions and
  deletions are O(n) at the splice site, not O(n²) string copies.
- Mutations are processed in ascending `cdna_start` order, with a running
  `cumulative_offset` carried forward. After a 5-bp insertion at c.100, the
  next mutation at c.200 is applied at index 100 + (200 − 1) + 5 = 304 of
  the sequence list, so coordinates remain correct on the running mutated
  sequence.
- **Reference-base sanity check.** When the MAF gives an explicit ref base,
  it must match the ORF base at that position; mismatches are reported with
  `status=skipped_ref_mismatch` (unless `--allow-ref-mismatch` is passed).
  When the ref was an `N` placeholder, this check is skipped — the ORF base
  is treated as the ground truth.
- Each variant type has its own splice rule (replace, delete, insert, or
  delins); the bookkeeping is in [mutations.py](mutations.py#L308-L460).

In the current pipeline, **each eligible variant is applied independently
to the wildtype ORF**, not stacked. That keeps the output variant-centred:
one FASTA entry tracks one mutation, and a peptide spanning the mutation
site is unambiguously evidence for that mutation.

---

## 9. Translation and effect classification

[`translate`](translate.py) reads the (possibly mutated) cDNA in frame from
index 0, translates triplet by triplet with the standard codon table, and
halts at the first stop codon by default. The reference ORF is also
translated once per `(sample, gene)` to provide a comparison protein.

Comparing the mutant protein to the reference protein, the effect is
classified in this priority order:

1. **`start_loss`** — first codon is no longer `ATG`.
2. **`frameshift`** — mutated ORF length differs from the reference by a
   non-multiple of 3.
3. **`stop_loss`** — reference protein had a stop and mutant translation
   ran past it (requires the ref stop position as input to disambiguate
   from an in-frame insertion).
4. **`nonsense`** — premature stop codon, mutant protein shorter than ref.
5. **`synonymous`** — protein sequences identical. *These are dropped from
   the FASTA* (logged as `sequence_status=skipped_synonymous` in
   `translation_report.tsv`). The MAF-level silent filter catches most;
   this catches the residual cases where the MAF said missense but the
   substitution turned out synonymous at the protein level.
6. **`in_frame_ins` / `in_frame_del`** — protein length changed by a
   multiple of 3 at the nucleotide level.
7. **`missense`** — same length, ≥1 amino-acid difference.

---

## 10. Extracting the variant-centred context

This is the step that makes the database **compact** rather than per-gene
full-length.

[`_variant_context_sequence`](build_mutated_protein_db.py) locates the
altered amino-acid interval by walking the reference and mutant protein
strings inward from both ends until they disagree
([`_altered_interval`](build_mutated_protein_db.py)). Then it pads the
interval by `--variant-context-flank-aa` (default 30 aa) on each side and
slices the mutant protein.

The flank of 30 aa is chosen so any HLA-I peptide of length 8–15 that
*overlaps* the altered residue is fully contained within the FASTA record.
For frameshifts and stop-losses, the novel downstream reading frame is
bounded at `--frameshift-tail-aa` (default 120 aa) to keep file size
manageable; in practice, biologically relevant peptides occur close to the
junction.

The FASTA header records the exact coordinates of the segment and the
altered sub-interval, e.g.:

```
>mutprot_3f9c1a2b8eef|TP53|ENST00000269305.9|mutated|sample=CSC0073.TXO|
 variants=c.743G>A|protein_change=p.R248Q|mutation_type=missense|
 coords=aa218-278;altered=aa248-249
```

The FASTA ID itself is a SHA-1 hash of `{variant_id}|{context_sequence}`.
Identical contexts therefore deduplicate naturally, and full traceability
(which sample(s), which variant(s), which transcript) is retained in the
sidecar `variant_to_protein_metadata.tsv`.

---

## 11. Peptide candidate generation (for downstream ML)

For each written context, [`_peptide_candidates`](build_mutated_protein_db.py)
slides a window of every length in `[--peptide-min-length,
--peptide-max-length]` (default 8–15) across the altered interval. Only
windows whose span *includes the altered residue* are emitted — peptides
from unchanged regions are not neoantigen-derived evidence and would just
duplicate reference proteome content. Windows containing `*` (stop) or
`X` (unknown) are dropped, and windows whose mutant sequence equals the
wildtype at the same coordinates are dropped.

If an HLA table is loaded, peptides are paired with each of the sample's
HLA alleles; otherwise one row is emitted with blank HLA. Binding
predictions and validated labels can be joined in via `--binding-tsv` and
`--labels-tsv`. The result is written as `ml_training_table.tsv`, one row
per `(peptide, HLA)` pair.

---

## 12. Outputs

After Step 01 completes, `mutation_reference/mutated_proteins/` contains:

| File | One row per | Purpose |
|---|---|---|
| `all_samples_mutated_proteins.fasta` | record | Cohort-level mutant search database |
| `per_sample/<SAMPLE>.mutated_proteins.fasta` | record | Sample-specific mutant search database |
| `mutation_application_report.tsv` | input variant | Did each variant apply? If not, why? |
| `translation_report.tsv` | mutated ORF | Lengths, frameshift flags, effect, FASTA ID |
| `variant_to_protein_metadata.tsv` | FASTA record | Full traceability from FASTA ID back to sample/gene/variant |
| `rna_feature_table.tsv` | variant | TPM, RNA VAF / depth / alt-count, filter status |
| `ml_training_table.tsv` | peptide × HLA | Mutant peptide candidates with WT counterpart and feature columns |
| `validation_report.tsv` | ORF warning | Issues detected in the reference ORF FASTA |

Step 03 ([03_merge_fragpipe_db.sh](03_merge_fragpipe_db.sh)) then unions
this mutant FASTA with the novel-transcript ORFs from Step 02 to produce
the final FragPipe input under `mutation_reference/final/`.

---

## 13. Common pitfalls and design notes

- **Variant-centred ≠ full mutant protein.** Peptides from unchanged
  regions of a mutated gene are *not* neoantigen evidence — they exist in
  the reference proteome already. Writing only the altered region prevents
  the search from being misled by trivial hits.
- **Independent variant application.** Compound effects of multiple
  mutations in the same ORF are not modelled in the current FASTA. If two
  variants in the same transcript would jointly produce a unique peptide
  that neither produces alone, that peptide is not currently in the search
  space. Adding a "stacked" mode is a tractable extension.
- **Longest-ORF tiebreak when transcript missing.** Variants without a
  matching transcript ID are mapped to the longest CDS for that gene. This
  is conservative for peptide identification but may not reflect the
  actually expressed isoform in the tumour. Transcript-specific behaviour
  could be improved with isoform-resolved RNA-seq later.
- **HLA pairing happens at the peptide table, not the FASTA.** The FASTA
  is HLA-agnostic — every sample's mutant contexts are in the search space
  regardless of HLA type. HLA restriction is applied downstream (peptide
  table for ML, binding prediction tools for ranking).
- **MS/MS coverage caveat.** The FASTA flank of 30 aa is sized for HLA-I
  peptides (8–15 aa). HLA-II peptides (commonly 13–25 aa) are still
  representable but the flank may need widening for very long peptides
  that include the altered residue near one terminus.
