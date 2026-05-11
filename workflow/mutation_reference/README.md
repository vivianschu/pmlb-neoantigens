# FragPipe Reference Database Pipeline

Builds compact sample-specific variant-centered mutant FASTA databases for FragPipe
non-specific HLA immunopeptidomics workflows, plus sidecar feature tables for
peptide/HLA machine-learning datasets.

For each sample, somatic mutations are mapped to cDNA ORF sequences at the
nucleotide level. Each applied coding variant is translated as a local mutant
protein context, and only the compact sequence around the altered amino-acid
interval is written to FASTA. Novel transcript proteins from RNA-seq can
optionally be merged in.
Silent MAF classifications and protein-identical translated products are
excluded from the mutant FASTA.

```
cDNA ORF FASTA
+ somatic MAF
                  →  map variant to ORF/codon  →  translate local context  →  compact mutant FASTA
+ novel transcript
  .pep (optional)
                  →  reformat headers            →             →  novel protein FASTA
                                                                          ↓
                                                             FragPipe database search
```

---

## Directory layout

```
workflow/mutation_reference/
├── 00_fetch_cdna_orfs.sh        SLURM: filter GENCODE → data/orf/cdna_orfs.fasta
├── 01_build_mutated_db.sh       SLURM: run build_mutated_protein_db.py
├── 02_novel_orfs.sh             SLURM array: TransDecoder ORFs → FragPipe FASTA
├── 03_merge_fragpipe_db.sh      SLURM: combine mutated + novel proteins
├── submit_build_db.sh           Top-level submission wrapper (chains all stages)
│
├── build_mutated_protein_db.py  Main CLI — builds variant-centered mutant contexts
├── orf.py                       Load and validate cDNA ORF FASTA records
├── mutations.py                 Parse MAF/HGVSc, apply mutations with offset tracking
├── translate.py                 Codon table, translation, mutation effect classification
├── fasta.py                     FragPipe-compatible FASTA header formatting and writing
│
└── tests/
    ├── test_apply_snv.py
    ├── test_apply_indel.py
    ├── test_frameshift_translation.py
    └── test_fasta_output.py
```

---

## Cluster paths

| Resource | Path |
|---|---|
| GENCODE v47 pc_transcripts | `/cluster/projects/livingbank/workspace/references/gencode.v47.pc_transcripts.fa.gz` |
| Novel transcript FASTAs | `/cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts/<SAMPLE>_novel_sequences.fasta` |
| TransDecoder .pep output | `/cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/novel_proteins/transdecoder_pep/` |
| TransDecoder work directories | `/cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/novel_proteins/transdecoder_work/` |
| Novel protein FASTAs (step 02 output) | `/cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/novel_proteins/` |
| Final merged database (step 03 output) | `/cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/final/` |

| Resource | Path |
|---|---|
| Somatic MAF | `/cluster/projects/livingbank/Project/Pan-organoid/Mutation/Org_exome.data_mutations_extended.gt4.202507.txt` |
| CNV | `/cluster/projects/livingbank/Project/Pan-organoid/CNV/data_CNV.txt` |
| cDNA ORF FASTA (step 00 output) | `data/orf/cdna_orfs.fasta` |
| Mutated proteins (step 01 output) | `/cluster/projects/livingbank/workspace/vivian/neo/mutation_reference/mutated_proteins/` |

---

## Running the pipeline

Run these scripts on the Slurm cluster. The Python environment used by the
jobs must already provide `pandas`; the job intentionally does not install
packages on compute nodes.

### Step 00 — Prepare the cDNA ORF reference (run once)

Filters the GENCODE v47 protein-coding transcript FASTA to ATG-starting ORFs and
reformats headers to `>GENE|TRANSCRIPT_ID|ORF`.

```bash
sbatch 00_fetch_cdna_orfs.sh
# Output: data/orf/cdna_orfs.fasta
```

### Step 01 — Build the variant-centered mutant database

Maps somatic mutations from the MAF file to each gene's cDNA ORF, translates the
single-variant mutant sequence, and writes per-sample and combined FASTAs
containing compact mutant protein contexts. These FASTA entries are not full
mutated genes or full transcript-derived proteins.

```bash
sbatch 01_build_mutated_db.sh
# Output: results/mutation_reference/mutated_proteins/
```

Submit from either the repository root or `workflow/mutation_reference/`. If
submitting from another directory, pass the source directory explicitly so Slurm
does not resolve helper scripts from `/var/spool/slurmd/`:

```bash
sbatch /path/to/workflow/mutation_reference/01_build_mutated_db.sh \
    --script-dir /path/to/workflow/mutation_reference
```

With a sample manifest to process a subset:

```bash
sbatch 01_build_mutated_db.sh \
    --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv
```

With optional RNA/HLA/binding inputs:

```bash
sbatch 01_build_mutated_db.sh \
    --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv \
    --expression-tsv /cluster/projects/livingbank/workspace/vivian/neo/rna/expression_tpm.tsv \
    --rna-variants-tsv /cluster/projects/livingbank/workspace/vivian/neo/rna/rna_variant_support.tsv \
    --binding-tsv /cluster/projects/livingbank/workspace/vivian/neo/binding/netmhcpan_predictions.tsv \
    --labels-tsv /cluster/projects/livingbank/workspace/vivian/neo/labels/peptide_labels.tsv \
    --binder-rank-threshold 2.0

# actual
sbatch 01_build_mutated_db.sh \
    --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv \
    --expression-tsv /cluster/projects/livingbank/Project/Pan-organoid/RNAseq/batch_corrected_TPM_all_genes.csv \
    --expression-gene-column gene_id \
    --hla-tsv /cluster/projects/livingbank/workspace/vivian/neo/hla/hla_summary.csv \
    --hla-sample-column sample
```

If present, `data/rna/hla_summary.csv` is used automatically for sample HLA
calls. You can override it with `--hla-tsv /path/to/hla.tsv`.

RNA filters are optional. If RNA data are missing and `--require-rna-filters`
is not set, the workflow keeps the variant and writes missing values in the
feature tables. To require RNA support:

```bash
sbatch 01_build_mutated_db.sh \
    --expression-tsv /path/to/expression_tpm.tsv \
    --rna-variants-tsv /path/to/rna_variant_support.tsv \
    --min-tpm 1 \
    --min-rna-vaf 0.05 \
    --min-rna-depth 10 \
    --require-rna-filters
```

### Step 02 — Predict and reformat novel transcript proteins (array job)

Runs TransDecoder for each per-sample novel transcript FASTA, caches the raw
TransDecoder `.pep` files under `novel_proteins/transdecoder_pep/`, keeps
TransDecoder work directories under `novel_proteins/transdecoder_work/`, and
writes per-sample FragPipe-ready `*.novel_proteins.fasta` files directly under
`novel_proteins/`.

```bash
N=$(find /cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts \
    -type f -name '*_novel_sequences.fasta' \
    ! -path '/cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts/logs/*' \
    | wc -l)
sbatch --array=0-$((N-1)) 02_novel_orfs.sh
# Output: neo/mutation_reference/novel_proteins/<SAMPLE>.novel_proteins.fasta
```

The array is indexed by actual `*_novel_sequences.fasta` files, not by
directories. Both flat files under `novel/novel_transcripts/` and nested
`<SAMPLE>/<SAMPLE>_novel_sequences.fasta` layouts are supported; sample IDs are
derived from the FASTA filename. Step 02 writes a per-sample status file under
`neo/mutation_reference/novel_proteins/status/`; samples with empty input or no
predicted ORFs get an empty `*.novel_proteins.fasta` plus a skipped status row.
The script submits Step 02 to the `himem` partition and requests 128G because
`TransDecoder.Predict` can exceed the older 16G or 64G allocation on large novel
FASTA inputs.

### Step 03 — Merge into final FragPipe database

Combines mutated protein FASTAs (step 01) and novel protein FASTAs (step 02) into
per-sample and combined databases.

```bash
sbatch 03_merge_fragpipe_db.sh
# Output: neo/mutation_reference/final/all_samples_combined.fasta
#         neo/mutation_reference/final/per_sample/<SAMPLE>.fragpipe_db.fasta
```

---

## Submitting the full pipeline in one command

`submit_build_db.sh` chains all stages with SLURM dependencies.

```bash
# Mutated proteins only (steps 00 must already be done)
bash submit_build_db.sh

# Full pipeline: mutated + novel + merge
bash submit_build_db.sh --with-novel

# Full pipeline with sample filter
bash submit_build_db.sh --with-novel \
    --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv

# Preview all sbatch commands without submitting
bash submit_build_db.sh --with-novel --dry-run

# Merge only (steps 01 and 02 already complete)
bash submit_build_db.sh --merge-only
```

Dependency graph:

```
[01] build_mutated_db  ──────────────────────┐
                                              ├─► [03] merge_fragpipe_db
[02] novel_orfs (array, runs in parallel) ───┘
```

---

## Output files

### Variant-centered mutant contexts (`results/mutation_reference/mutated_proteins/`)

| File | Description |
|---|---|
| `all_samples_mutated_proteins.fasta` | Compact variant-centered mutant contexts, all samples combined |
| `per_sample/<SAMPLE>.mutated_proteins.fasta` | Compact variant-centered mutant contexts for one sample |
| `mutation_application_report.tsv` | One row per mutation — applied or skipped with reason |
| `translation_report.tsv` | One row per mutated ORF — lengths, frameshift, effect type |
| `validation_report.tsv` | ORF validation warnings (missing ATG, length, internal stops) |
| `variant_to_protein_metadata.tsv` | FASTA ID to sample/gene/transcript/variant/protein-effect sidecar table |
| `rna_feature_table.tsv` | Per-variant TPM/RNA VAF/RNA depth/RNA alt-read features and filter status |
| `ml_training_table.tsv` | Peptide-level mutant candidates for ML training |

Use `--compress-tables` to write these TSV outputs as `.tsv.gz`.

The mutant FASTA is intentionally compact. A mutation-associated peptide must
directly contain the altered amino acid sequence, span the altered indel
junction, or originate from the novel downstream reading frame created by a
frameshift. Peptides from an unchanged region of a mutated gene are not treated
as mutation-derived evidence.

### Final database (`neo/mutation_reference/final/`)

| File | Description |
|---|---|
| `all_samples_combined.fasta` | Mutated + novel proteins, all samples |
| `per_sample/<SAMPLE>.fragpipe_db.fasta` | Per-sample FragPipe input |
| `merge_summary.tsv` | Per-sample record counts (mutated / novel / total) |

### FASTA header format

```
>mutprot_<sha1>|GENE|TRANSCRIPT_ID|mutated|sample=SAMPLE_ID|variants=c.743G>A|protein_change=p.R248Q|mutation_type=missense|coords=aa218-278;altered=aa248-249
>SAMPLE_ID|GENE|TRANSCRIPT_ID|reference                      (--include-reference only)
>SAMPLE_ID|TRANSCRIPT_ID|novel|source=transdecoder|orf_type=complete|strand=+
```

Identical compact mutant context sequences are deduplicated by default. The
FASTA ID is sequence-derived, and all sample/variant traceability is retained in
`variant_to_protein_metadata.tsv`. Disable this only for debugging with
`--no-deduplicate`.

## Optional RNA/HLA/input table formats

Default column names can be overridden on the Python CLI.

Expression table:

| Column | Meaning |
|---|---|
| `sample_id` | Sample barcode matching the MAF after manifest normalization |
| `gene` | Gene symbol |
| `transcript_id` | Transcript ID; may be empty if gene-level TPM is used |
| `TPM` | Gene or transcript TPM |

Long expression tables may be TSV or CSV. Wide gene-by-sample matrices are also
accepted; for example `batch_corrected_TPM_all_genes.csv` uses `gene_id` as the
first column and sample IDs as the remaining columns:

```bash
--expression-tsv /cluster/projects/livingbank/Project/Pan-organoid/RNAseq/batch_corrected_TPM_all_genes.csv \
--expression-gene-column gene_id
```

Sample suffixes such as `_pmlb` and `_novo` are normalized. If multiple RNA
columns collapse to the same canonical sample ID, the maximum TPM is used for
that sample/gene. The expression gene key must match the mutation/ORF gene key
for TPM to populate; if mutations use Hugo symbols and expression uses Ensembl
IDs only, a symbol-keyed expression table or mapping step is needed.

RNA variant support table:

| Column | Meaning |
|---|---|
| `sample_id` | Sample barcode |
| `gene` | Gene symbol |
| `transcript_id` | Optional transcript ID |
| `variant_id` or `HGVSc` | Variant key matching MAF-derived IDs or HGVSc |
| `rna_vaf` | RNA variant allele fraction |
| `rna_depth` | RNA read depth |
| `rna_alt_count` | RNA alternate-read count |

HLA table:

| Column | Meaning |
|---|---|
| `sample_id` | Sample barcode |
| `hla_allele` | One allele or comma/semicolon/space-separated alleles |

The workflow also accepts the repository RNA HLA summary CSV directly:
`data/rna/hla_summary.csv` with columns `sample,status,reason,hla_a1,hla_a2,
hla_b1,hla_b2,hla_c1,hla_c2`. Calls with non-success status are skipped,
sample suffixes such as `_novo` and `_pmlb` are normalized, duplicate alleles
are collapsed, and alleles are reported in two-field `HLA-A*02:01` style.

Binding and label tables are optional joins keyed by sample, mutant peptide, and
HLA allele. Labels are never fabricated. `presented_label` should only be filled
from appropriate peptide evidence such as accepted FragPipe PSM/peptide evidence,
and `immunogenic_label` should only come from validated assays. If
`--binder-ic50-threshold` or `--binder-rank-threshold` is supplied, missing
`binder_label` values are derived from that explicit threshold.

## ML training table

`ml_training_table.tsv` has one row per mutant peptide candidate and HLA allele
when HLA alleles are available. If HLA input is absent, one row is emitted with
blank HLA fields. Peptides are generated only from windows overlapping the
altered amino-acid interval, for lengths controlled by `--peptide-min-length`
and `--peptide-max-length` (defaults 8-15). Long frameshift/stop-loss tails are
bounded by `--frameshift-tail-aa` to limit file size. The FASTA context around
each altered interval is controlled by `--variant-context-flank-aa` (default 30
aa on each side).

Core columns include mutant peptide, corresponding wildtype peptide where
available, HLA allele/key, TPM, RNA VAF/depth/alt-count, DNA VAF/depth/alt-count,
binding IC50/rank, peptide length, gene, transcript/protein ID, variant ID,
amino-acid change, mutation type, sample ID, and nullable binder/presented/
immunogenic labels.

## FragPipe HLA search assumptions

Use the per-sample FASTA from `final/per_sample/<SAMPLE>.fragpipe_db.fasta` for
sample-specific searches whenever possible. Use the cohort FASTA only when a
cohort-level search space is intentional. Recommended FragPipe assumptions for
HLA immunopeptidomics are non-specific or unspecific digestion, no enzyme
cleavage constraints, peptide lengths matching the HLA class under study
commonly 8-15 for class I, and no artificial requirement for tryptic termini.
Keep the target-decoy and contaminant handling consistent with the rest of the
project's FragPipe configuration.

---

## Direct CLI validation

```bash
python workflow/mutation_reference/build_mutated_protein_db.py \
    --orf-fasta  data/orf/cdna_orfs.fasta \
    --mutations  /cluster/projects/livingbank/Project/Pan-organoid/Mutation/Org_exome.data_mutations_extended.gt4.202507.txt \
    --outdir     results/mutation_reference/mutated_proteins \
    --per-sample
```

Key options for the Python CLI:

| Flag | Default | Description |
|---|---|---|
| `--per-sample` | off | Write one FASTA per sample |
| `--include-reference` | off | Also write unmodified reference proteins |
| `--allow-ref-mismatch` | off | Apply mutations even when HGVSc ref ≠ ORF base |
| `--translate-through-stop` | off | Continue translation past stop codons |
| `--sample-manifest FILE` | — | Process only samples in FILE; accepts a plain list or a TSV with a `sample` column (e.g. `260313_manifest_final.tsv`); `_pmlb`/`_novo` suffixes are stripped automatically to match MAF barcodes |
| `--expression-tsv FILE` | — | Optional TPM table |
| `--rna-variants-tsv FILE` | — | Optional RNA VAF/depth/alt-read table |
| `--hla-tsv FILE` | — | Optional sample HLA table |
| `--binding-tsv FILE` | — | Optional binding predictions |
| `--labels-tsv FILE` | — | Optional validated labels or peptide-evidence labels |
| `--compress-tables` | off | Write sidecar/report TSVs as `.tsv.gz` |
| `--variant-context-flank-aa N` | 30 | Flanking amino acids to include around each altered interval in mutant FASTA records |
| `--frameshift-tail-aa N` | 120 | Maximum novel downstream frameshift/stop-loss amino-acid region considered for FASTA contexts and peptide candidates |

---

## Mutation handling

Mutations are parsed from HGVSc notation and applied to the cDNA ORF at the
nucleotide level. The MAF file uses `N` as a reference-base placeholder
(e.g. `c.157delNNNNN`, `c.10138N>A`); the pipeline reads the actual reference
from the ORF sequence at that position.

| Variant type | Example HGVSc | Behaviour |
|---|---|---|
| SNV | `c.743G>A` | Replace single base |
| MNV | `c.743_744GC>TT` | Replace multiple bases |
| DEL | `c.157_161delNNNNN` | Delete bases; length from coordinate span |
| INS | `c.743_744insAA` | Insert after position 743 |
| INDEL | `c.743_745delinsAA` | Replace ref bases with alt bases |

Mutations are applied in ascending coordinate order. Insertions and deletions
shift all subsequent mutation coordinates automatically.

Translated effects reported: `missense`, `synonymous`, `nonsense`, `frameshift`,
`in_frame_ins`, `in_frame_del`, `stop_loss`, `start_loss`.

Synonymous translated effects are reported as `skipped_synonymous` in
`translation_report.tsv` and are not written to the mutant FASTA or peptide
candidate table.

For FASTA generation, each eligible variant is applied and translated
independently so the output remains variant-centered. Missense records retain
only contexts containing the mutated residue, in-frame insertion/deletion
records retain contexts containing the altered stretch or novel junction, and
frameshift records retain the novel downstream reading-frame context plus only
the configured upstream flank needed for peptide digestion/search. The full
unchanged upstream protein is not written as mutant evidence.

---

## Tests

```bash
cd workflow/mutation_reference
python tests/test_apply_snv.py
python tests/test_apply_indel.py
python tests/test_frameshift_translation.py
python tests/test_fasta_output.py
```

Requires `pandas` (`pip install pandas`).
