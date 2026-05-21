# Variant Calling and Benchmarking of HG002 Whole-Exome Sequencing Data

A reproducible variant calling and benchmarking pipeline applied to the GIAB HG002 whole-exome sequencing sample on chromosome 22, with stratified evaluation across four coverage levels (2x, 10x, 40x, 80x). Completed as a group project for the Postgraduate Bioinformatics Diploma at Nile University.

**Authors:** Mohammed A. Ahmed, Yasmine M. Elgohary, Ahmad K. Ahmad

## Question

How does sequencing depth affect variant calling accuracy in whole-exome sequencing data, and at what coverage level do gains in precision, recall, and genotype concordance start to saturate? Answering this on a well-characterized truth sample (GIAB HG002) with properly exome-restricted depth estimation provides a practical reference for cost-effective study design.

## Pipeline

The pipeline starts from the GIAB HG002 WES BAM file (originally aligned to an earlier reference build) and produces stratified benchmarking metrics on chr22. The BAM is first name-sorted to preserve read pairing during BAM-to-FASTQ conversion (a critical step identified after observing read-pair loss with direct conversion). FASTQs are then realigned to GRCh38 with BWA-MEM, deduplicated using Picard MarkDuplicates, and subset to chr22. Coverage is estimated using the Agilent SureSelect V5 exome capture BED file: a key methodological point is that genome-wide depth on chr22 (~14x) substantially underrepresents effective WES coverage, while exome-restricted depth (~215x across ~1.09 Mb of targeted bases) reflects the true sequencing depth and motivates a realistic downsampling strategy. The BAM is downsampled to ~80x, ~40x, ~10x, and ~2x using Picard DownsampleSam with fixed random seeds. Each coverage level undergoes Base Quality Score Recalibration (GATK BaseRecalibrator and ApplyBQSR with dbSNP138 as known sites) and variant calling with GATK HaplotypeCaller in GVCF mode, followed by joint genotyping with GenotypeGVCFs. Variants are stratified into SNPs and INDELs and into known (rsID-annotated) versus novel subsets, then assessed via density plots across the standard GATK quality metrics (QD, MQ, MQRankSum, FS, SOR, ReadPosRankSum, DP, InbreedingCoeff). Coverage-aware hard filters are applied with GATK VariantFiltration, including a per-coverage maximum DP threshold computed as mean(DP) + 5*SD(DP) to account for depth-dependent scaling. PASS-only VCFs are benchmarked against the GIAB HG002 GRCh38 truth set restricted to chr22 and intersected with the exome capture BED, using hap.py with the vcfeval engine in a dedicated Conda environment.

## Results

Benchmarking shows the expected accuracy-versus-depth curve. SNP F1 score rises sharply from 0.31 at 2x to 0.80 at 10x and 0.94 at 40x, then plateaus by 80x. INDELs lag SNPs at every coverage level (F1 of 0.21, 0.62, 0.77, 0.87 across the four depths). False positives per Mb drop steeply with depth (104, 38, 11, 10 for SNPs; 5, 17, 11, 6 for INDELs). Genotype concordance is high once coverage exceeds 10x, reaching 0.9988 for SNPs and 0.9833 for INDELs at 80x. The combined picture shows that the marginal benefit of sequencing beyond ~40x is small for chr22 WES on this sample, which is directly useful when budgeting for WES studies. Density-plot exploration also showed that novel (non-rsID) variants are enriched in the low-quality tails of every metric, validating coverage-aware hard filtering as a meaningful preprocessing step.

## Tools

BWA-MEM for alignment; Picard (MarkDuplicates, DownsampleSam); samtools; GATK4 (BaseRecalibrator, ApplyBQSR, AnalyzeCovariates, HaplotypeCaller, GenotypeGVCFs, SelectVariants, VariantFiltration); bedtools for region intersection; hap.py with the RTG vcfeval engine for benchmarking; R for visualization. Reference data: GRCh38, dbSNP138, GIAB HG002 truth VCF and confident regions BED, and Agilent SureSelect V5 capture BED.

## Repository contents

- `Project_Script_Group5.sh`: end-to-end bash pipeline from alignment through benchmarking
- `benchmarking_summary_chr22.tsv`: unified benchmarking metrics across all coverages and variant types
- `Density_plots.zip`: variant quality metric density plots stratified by coverage, type, and novelty
- `SNP_INDEL_filter_reports.zip`: per-coverage filtering summaries
- `benchmarking_plots.zip`: precision, recall, F1, genotype concordance, and FP per Mb plots from hap.py
- `HG002_WES_Variant_Calling_and_Benchmarking_Workflow.pdf`: visual workflow infographic
- `Variant_Calling_And_Benchmarking_Of_Hg002_WES_Data.pdf`: full scientific report
- `Technical_Journey_-_Variant_Calling_and_Benchmarking_Pipeline_Development.pdf`: technical log of decisions, corrections, and iterations

## Data

Source data and references:
- GIAB HG002 WES BAM (input sample)
- GIAB HG002 GRCh38 truth set and high-confidence regions (benchmarking)
- GRCh38 reference genome
- dbSNP138
- Agilent SureSelect V5 exome capture BED

None of these are redistributed in this repository. Source URLs are documented in the pipeline script.

## Reproducibility

The pipeline assumes a Linux environment with Conda. A primary Conda environment is used for the main pipeline (samtools, BWA, Picard, GATK4, bedtools); a separate environment (`happyenv`) is required for hap.py because it depends on Python 2.7. Both environment setups are documented in `Project_Script_Group5.sh`. Running the full workflow on chr22 from raw BAM to benchmarked PASS-only VCFs takes a few hours on a standard workstation.
