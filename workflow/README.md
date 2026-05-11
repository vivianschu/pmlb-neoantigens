# PMLB Neoantigen Discovery Pipeline

Multi-cancer organoid neoantigen discovery pipeline. Processes HLA typing, novel transcript
extraction, somatic mutation application, and SNV/indel filtering into a unified, ranked
neoantigen candidate list per sample.

This pipeline is designed for **tumor-only** samples — no matched normal tissue is available.

---

## Pipeline overview

```
RNA BAMs  ──► 01 HLA typing (arcasHLA)           ──► hla_summary.csv
          ──► 02 Novel transcript extraction       ──► per-sample FASTA / GTF
          ──► 03 TransDecoder ORF prediction       ──► per-sample .transdecoder.pep

cDNA ORF FASTA + somatic MAF
          ──► mutation_reference/01 Build mutated DB   ──► mutated_proteins/
          ──► mutation_reference/02 Novel ORFs         ──► novel_proteins/
          ──► mutation_reference/03 Merge              ──► final FragPipe DB

MAF + CNV + TPM + HLA
          ──► snv/   SNV/indel filtering          ──► results/snv/      [planned]
          ──► cnv/   CNV cohort analysis           ──► results/cnv/      [planned]
          ──► novel/ Novel peptide tiling          ──► results/novel/    [planned]
          ──► integrate/ Unified integration       ──► results/integrated/ [planned]
```

Steps 01–03 (RNA) and the mutation_reference pipeline run on the **cluster** via SLURM.
SNV/CNV/novel/integration steps run **locally**.

---

## Data layout

```
data/
├── mutation/
│   └── Org_exome.data_mutations_extended.gt4.202507.txt  # somatic MAF (VEP-annotated)
├── cnv/
│   └── data_CNV.txt                                       # Hugo_Symbol × samples, -2..2
├── rna/
│   ├── batch_corrected_TPM_all_genes.csv                  # Ensembl gene IDs × samples
│   ├── hla_summary.csv                                    # arcasHLA typed alleles
│   ├── metadata.csv                                       # sample → cancer_type, batch
│   └── 260313_manifest_final.tsv                          # SLURM job manifest
├── reference/                                             # reference genome/annotation files
└── atac/                                                  # ATAC-seq data
```

Cluster source files:
```
/cluster/projects/livingbank/Project/Pan-organoid/
├── Mutation/Org_exome.data_mutations_extended.gt4.202507.txt
└── CNV/data_CNV.txt
```

---

## Cluster directory structure

```
/cluster/projects/livingbank/workspace/vivian/neo/
├── 260313_manifest_final.tsv
├── hla/                                 # Step 01 outputs
│   ├── logs/
│   └── <sample>/
│       ├── <sample>.genotype.json
│       └── <sample>.genotype.tsv
├── novel/                               # Step 02 outputs
│   ├── logs/
│   ├── novel_translated/                # TransDecoder .pep files
│   ├── all_novel_sequences.fasta
│   ├── all_novel_transcripts.gtf
│   ├── novel_transcripts_summary.csv
│   └── <sample>/
│       ├── <sample>_gffcmp.annotated.gtf
│       ├── <sample>_novel_transcripts.gtf
│       └── <sample>_novel_sequences.fasta
└── mutation_reference/                  # mutation_reference pipeline outputs
    ├── logs/
    ├── orf/
    │   └── cdna_orfs.fasta              # Step 00 output
    ├── mutated_proteins/                # Step 01 output
    │   ├── all_samples_mutated_proteins.fasta
    │   ├── per_sample/<SAMPLE>.mutated_proteins.fasta
    │   ├── mutation_application_report.tsv
    │   ├── translation_report.tsv
    │   └── validation_report.tsv
    ├── novel_proteins/                  # Step 02 output
    │   └── <SAMPLE>.novel_proteins.fasta
    └── final/                           # Step 03 output (FragPipe DB)
        ├── all_samples_combined.fasta
        ├── per_sample/<SAMPLE>.fragpipe_db.fasta
        └── merge_summary.tsv
```

Input BAMs (STAR-aligned) and StringTie GTFs:
```
/cluster/projects/livingbank/workspace/vivian/RNA/process/rna/<sample>/
    <sample>_Aligned.sortedByCoord.out.bam
    <sample>_annotation.gtf    (StringTie output, required for Step 02)
```

---

## Cluster prerequisites

- SLURM account: `hansengroup`
- Modules: `kallisto`, `samtools`, `gffread`, `gffcompare`
- Tools in PATH: `arcasHLA`, `stringtie2`
  - stringtie2 / gffcompare: `/cluster/home/t117036uhn/local_soft/`
- Reference: `hg38_ek12` at `/cluster/projects/livingbank/workspace/references/hg38_ek12/`
- Reference GTF: `/cluster/projects/livingbank/workspace/references/hg_38/gencode.v25.annotation.gtf`
- GENCODE v47 pc_transcripts: `/cluster/projects/livingbank/workspace/references/gencode.v47.pc_transcripts.fa.gz`

---

## Step 01 — HLA Typing (cluster) · `rna/01_hla_typing/`

Runs [arcasHLA](https://github.com/RabadanLab/arcasHLA) on each sample's STAR-aligned BAM
to type HLA-A, -B, -C alleles at two-field resolution.

**Fallback:** If arcasHLA fails, the sample is assigned population-frequency priors
(`HLA-A02:01`, `HLA-A24:02`, `HLA-B07:02`, `HLA-B44:02`, `HLA-C07:01`, `HLA-C07:02`).
Imputed samples are flagged with `hla_imputed=True`.

### Scripts

| Script | Purpose |
|---|---|
| `hla_typing_array.sh` | SLURM array — one task per sample |
| `submit_hla.sh` | Submission wrapper: validates inputs, submits array, merges results |

### Run

```bash
conda activate neo
./submit_hla.sh --dry-run    # check inputs
./submit_hla.sh              # submit jobs
./submit_hla.sh --merge-only # merge after completion
```

### Outputs

| File | Description |
|---|---|
| `<sample>/<sample>.genotype.json` | arcasHLA raw genotype |
| `<sample>/<sample>.genotype.tsv` | TSV version |
| `hla_summary.csv` | One row per sample: HLA-A/B/C alleles and status |
| `failed_samples.txt` | Samples without a genotype, with reason |

### Monitor

```bash
squeue -j <JOB_ID>
tail -f /cluster/projects/livingbank/workspace/vivian/neo/hla/logs/hla_<JOBID>_<TASKID>.out
```

---

## Step 02 — Novel Transcript Extraction (cluster) · `rna/02_extract_novel/`

Annotates per-sample StringTie assemblies against GENCODE v25 with gffcompare, extracts
novel transcripts (class codes `u`, `n`, `j`), and exports sequences with gffread.

| Class code | Description |
|---|---|
| `u` | Intergenic — novel locus, no overlap with any reference gene |
| `n` | Intronic — transcript fully within a known gene intron |
| `j` | Novel junction — known gene with at least one novel splice site |

### Scripts

| Script | Purpose |
|---|---|
| `extract_novel.sh` | SLURM array — one task per sample |
| `submit_extract_novel.sh` | Submission wrapper |

### Run

```bash
./submit_extract_novel.sh --dry-run
./submit_extract_novel.sh
./submit_extract_novel.sh --merge-only
```

### Outputs

Per sample (`novel/<sample>/`): annotated GTF, filtered GTF, novel FASTA.

Collected (`novel/`): `all_novel_sequences.fasta`, `all_novel_transcripts.gtf`,
`novel_transcripts_summary.csv`.

---

## Step 03 — TransDecoder ORF Prediction (cluster) · `rna/03_translate/`

Runs [TransDecoder](https://github.com/TransDecoder/TransDecoder) on the per-sample novel
FASTAs to produce `.transdecoder.pep` files. Minimum ORF length is 25 aa (TransDecoder
default is 100 aa, which is too strict for short junction transcripts).

### Run

```bash
N=$(find /cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts \
    -type f -name '*_novel_sequences.fasta' \
    ! -path '/cluster/projects/livingbank/workspace/vivian/neo/novel/novel_transcripts/logs/*' \
    | wc -l)

sbatch --array=0-$((N-1)) workflow/rna/03_translate/run_transdecoder.sh

# Check one array task without running TransDecoder
SLURM_ARRAY_TASK_ID=0 bash workflow/rna/03_translate/run_transdecoder.sh --dry-run
```

The array is indexed by actual `*_novel_sequences.fasta` files and supports both
flat files under `novel_transcripts/` and nested
`novel_transcripts/<sample>/<sample>_novel_sequences.fasta` layouts. The job
requests 64G because `TransDecoder.Predict` can exceed smaller allocations on
large per-sample novel FASTAs.

### Outputs

Written alongside each FASTA in `novel_translated/`:

| File | Description |
|---|---|
| `<sample>_novel_sequences.fasta.transdecoder.pep` | Predicted protein sequences |
| `<sample>_transdecoder_wd/` | TransDecoder working directory |

---

## mutation_reference pipeline (cluster) · `mutation_reference/`

Builds sample-specific variant-centered mutant and novel protein FASTA databases for FragPipe
immunopeptidomics searches. See `mutation_reference/README.md` for full documentation.

### Overview

```
cDNA ORF FASTA + somatic MAF  →  map variants to ORF/codon → translate local context → compact mutant FASTA
novel transcript .pep          →  reformat headers            →              →  novel protein FASTA
                                                                                        ↓
                                                                           FragPipe database search
```

### Scripts

| Script | Purpose |
|---|---|
| `00_fetch_cdna_orfs.sh` | Extract CDS sequences from GENCODE v47 → `cdna_orfs.fasta` |
| `01_build_mutated_db.sh` | Map somatic mutations, translate local altered contexts, write compact mutant FASTAs |
| `02_novel_orfs.sh` | SLURM array: reformat TransDecoder .pep → FragPipe FASTA |
| `03_merge_fragpipe_db.sh` | Combine mutated + novel into final per-sample DB |
| `submit_build_db.sh` | Top-level wrapper: chains all stages with SLURM dependencies |

### Run

```bash
cd /cluster/projects/livingbank/workspace/vivian/neo/mutation_reference

# Dry run to verify paths
bash submit_build_db.sh --with-novel --dry-run

# Full pipeline: mutated + novel + merge
bash submit_build_db.sh --with-novel \
    --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv

# Mutated proteins only (no novel)
bash submit_build_db.sh \
    --sample-manifest /cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv

# Merge only (steps 01 and 02 already done)
bash submit_build_db.sh --merge-only
```

### Sample manifest note

`--sample-manifest` accepts `260313_manifest_final.tsv` directly. The `_pmlb`/`_novo`
suffixes in the `sample` column are stripped automatically to match MAF
`Tumor_Sample_Barcode` values.

---

## Design rationale: tumor-only sequencing

Without a matched normal sample, germline contamination problems arise that don't exist in
paired tumor/normal analysis. The pipeline addresses these through a **somatic confidence
score** incorporating VAF, depth, gnomAD frequency, gene-level germline risk, and RNA
support.

**Three-axis scoring model** (planned for SNV/integration steps):

```
somatic_confidence  ×  bio_opportunity  ×  mhc_presentation
     (0..1)             (0..∞, pre-MHC)      (0..1, post-MHC)
```

- **somatic_confidence**: how likely is this variant truly somatic?
- **bio_opportunity**: how much mutant peptide substrate is available for MHC loading?
- **mhc_presentation**: how well does the mutant peptide bind and distinguish from wildtype?

---

## Sample naming conventions

Sample IDs use a base format (e.g., `PPTO0002.TPO`) with optional batch/replicate suffixes:

| Suffix | Meaning |
|---|---|
| `_pmlb` | Primary PMLB batch |
| `_novo`, `_novo1`, `_novo2` | Novogene sequencing batches |
| `.r1`, `.r2` | Replicate sequencing runs |
| `.NPT` / `.NPO` | Normal patient tissue / organoid (excluded from tumor analysis) |

Sample type codes:

| Code | Type |
|---|---|
| `.TPO` | Tumor patient organoid |
| `.TPT` | Tumor patient tissue |
| `.TXO` | Tumor xenograft organoid |
| `.NPT` | Normal patient tissue |
| `.NPO` | Normal patient organoid |

All scripts normalize batch suffixes: `PPTO0002.TPO_pmlb` → `PPTO0002.TPO` before
cross-file joining.

---

## Sample manifest

Both cluster steps derive the sample list from:
```
/cluster/projects/livingbank/workspace/vivian/neo/260313_manifest_final.tsv
```
Columns: `sample`, `read1`, `read2`, `source_sample_names`.
To add or remove samples, edit the manifest — no script changes needed.
