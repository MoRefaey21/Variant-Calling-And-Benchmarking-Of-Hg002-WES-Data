#!/usr/bin/env bash

# -------------------------
# Setting environment (paths and directories)
# -------------------------

# Activating the conda environment
conda activate ngs1

# Project root path
PROJ_DIR="$HOME/workdir/GATK_tutorial/project"

# Creating the project directory tree
mkdir -p "$PROJ_DIR"

# Creating standard subdirectories used later
mkdir -p \
  "$PROJ_DIR/filters" \
  "$PROJ_DIR/benchmarking/truth" \
  "$PROJ_DIR/benchmarking/bed" \
  "$PROJ_DIR/benchmarking/tmp" \
  "$PROJ_DIR/benchmarking/results"

# Moving to the project root directory
cd "$PROJ_DIR" || exit 1

# -------------------------
# Downloading the WES BAM file
# -------------------------

wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

# -------------------------
# Converting BAM to FASTQ (for re-alignment to GRCh38)
# -------------------------

# Name-sorting the BAM so paired reads stay together during BAM to FASTQ conversion
samtools sort -n -@14 \
  -o HG002_namesorted.bam \
  151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

# Converting to paired FASTQs (discarding unpaired outputs)
samtools fastq -@14 \
  -1 HG002_R1.fastq.gz \
  -2 HG002_R2.fastq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n \
  HG002_namesorted.bam

# -------------------------
# Downloading and indexing the reference genome (GRCh38)
# -------------------------

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

# Samtools FASTA index
samtools faidx Homo_sapiens_assembly38.fasta

# Picard sequence dictionary
java -Xmx30g -jar $picard_path/picard.jar CreateSequenceDictionary \
    R=Homo_sapiens_assembly38.fasta \
    O=Homo_sapiens_assembly38.dict

# BWA index
bwa index Homo_sapiens_assembly38.fasta

# -------------------------
# Aligning FASTQs to GRCh38
# -------------------------

REF="Homo_sapiens_assembly38.fasta"
R1="HG002_R1.fastq.gz"
R2="HG002_R2.fastq.gz"

bwa mem -t 14 -M \
  -R "@RG\tID:HG002_GRCh38\tSM:HG002\tPL:ILLUMINA\tLB:HG002_WES" \
  "$REF" "$R1" "$R2" \
  > HG002_GRCh38.sam

# Converting SAM to BAM and sorting
samtools view -@14 -b HG002_GRCh38.sam \
  | samtools sort -@14 -o HG002_GRCh38.sorted.bam

# Indexing the BAM file
samtools index -@14 HG002_GRCh38.sorted.bam

# -------------------------
# Marking duplicates
# -------------------------

picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0
java -Xmx16g -jar $picard_path/picard.jar MarkDuplicates \
  INPUT=HG002_GRCh38.sorted.bam \
  OUTPUT=HG002_GRCh38.dedup.bam \
  METRICS_FILE=HG002_GRCh38.metrics.txt \
  CREATE_INDEX=true

# -------------------------
# Subsetting chr22 from the deduplicated BAM
# -------------------------

BAM_DEDUP="HG002_GRCh38.dedup.bam"

samtools view -@14 -b "$BAM_DEDUP" chr22 > HG002_22.bam
samtools sort -@14 -o HG002_22.sorted.bam HG002_22.bam
samtools index -@14 HG002_22.sorted.bam

# Checking the depth

samtools depth -@14 -b S04380110_Regions.sorted.bed HG002_22.sorted.bam \
  | awk '{sum+=$3} END {print "Average =",sum/NR,"No of covered Nuc =",NR}'

# Average = 215.42 No of covered Nuc = 1092045

# -------------------------
# Subsampling to 80x, 40x, 10x, and 2x using Picard DownsampleSam
# -------------------------

# 80×: P80 = 80/215.42 ≈ 0.3713
# 40×: P40 = 40/215.42 ≈ 0.1856
# 10×: P10 = 10/215.42 ≈ 0.0464
# 2×:  P2  = 2/215.42  ≈ 0.0093

picard_path=$CONDA_PREFIX/share/picard-* ## 2.21.7-0

# 80X
java -Xmx30g -jar $picard_path/picard.jar DownsampleSam \
  I=HG002_22.sorted.bam \
  O=HG002_22_80x.bam \
  P=0.3713 \
  RANDOM_SEED=42 \
  STRATEGY=ConstantMemory

# 40X
java -Xmx30g -jar $picard_path/picard.jar DownsampleSam \
  I=HG002_22.sorted.bam \
  O=HG002_22_40x.bam \
  P=0.1856 \
  RANDOM_SEED=43 \
  STRATEGY=ConstantMemory

# 10X
java -Xmx30g -jar $picard_path/picard.jar DownsampleSam \
  I=HG002_22.sorted.bam \
  O=HG002_22_10x.bam \
  P=0.0464 \
  RANDOM_SEED=44 \
  STRATEGY=ConstantMemory

# 2X
java -Xmx30g -jar $picard_path/picard.jar DownsampleSam \
  I=HG002_22.sorted.bam \
  O=HG002_22_2x.bam \
  P=0.0093 \
  RANDOM_SEED=45 \
  STRATEGY=ConstantMemory

# Sorting the downsampled BAMs
for bam in *x.bam
do
  samtools sort -@14 "$bam" -o "${bam%.bam}.sorted.bam"
done

# Indexing the downsampled BAMs
for bam in *x.sorted.bam
do
  samtools index -@ 14 "$bam"
done

# Checking the depth for each downsampled BAM
for bam in *x.sorted.bam; do 
    samtools depth -a -b S04380110_Regions.sorted.bed "$bam" \
    | awk -v fname="$bam" '{sum+=$3} END {print fname, "Average =", sum/NR, "No of covered Nuc =", NR}'
done

HG002_22_10x.sorted.bam Average = 9.99063 No of covered Nuc = 1092705
HG002_22_2x.sorted.bam Average = 2.02348 No of covered Nuc = 1092705
HG002_22_40x.sorted.bam Average = 39.924 No of covered Nuc = 1092705
HG002_22_80x.sorted.bam Average = 79.8469 No of covered Nuc = 1092705

# -------------------------
# Downloading known variants
# -------------------------

wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

# Indexing
gatk IndexFeatureFile -I Homo_sapiens_assembly38.dbsnp138.vcf

# -------------------------
# BQSR
# -------------------------

for sample in *.sorted.bam; do
    name=${sample%.sorted.bam}

    gatk --java-options "-Xmx30G" BaseRecalibrator \
        -R Homo_sapiens_assembly38.fasta \
        -I $sample \
        --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
        -O $name.report

    gatk --java-options "-Xmx30G" ApplyBQSR \
        -R Homo_sapiens_assembly38.fasta \
        -I $sample \
        -bqsr $name.report \
        -O $name.bqsr.bam

    gatk --java-options "-Xmx30G" BaseRecalibrator \
        -R Homo_sapiens_assembly38.fasta \
        -I $name.bqsr.bam \
        --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
        -O $name.report2

    gatk --java-options "-Xmx30G" AnalyzeCovariates \
        -before $name.report \
        -after $name.report2 \
        -plots $name.pdf
done

# -------------------------
# HaplotypeCaller
# -------------------------

for sample in *.bqsr.bam; do
    name=${sample%.bqsr.bam}

    gatk --java-options "-Xmx30G" HaplotypeCaller \
        -R Homo_sapiens_assembly38.fasta \
        -I $sample \
        -L chr22 \
        --emit-ref-confidence GVCF \
        --pcr-indel-model NONE \
        --native-pair-hmm-threads 14 \
        -O $name.gvcf
done

# -------------------------
# Genotyping all GVCFs
# -------------------------

for gvcf in *.gvcf; do
    name=${gvcf%.gvcf}

    gatk --java-options "-Xmx30G" GenotypeGVCFs \
        -R Homo_sapiens_assembly38.fasta \
        -V "$gvcf" \
        -L chr22 \
        --max-alternate-alleles 6 \
        --dbsnp Homo_sapiens_assembly38.dbsnp138.vcf \
        -O "${name}_raw_variants_ann.vcf"
done

# -------------------------
# Checking how many variants got annotated
# -------------------------

for vcf in *_raw_variants_ann.vcf; do
  echo "==> $vcf"
  grep -v "^#" "$vcf" | awk '{print $3}' | grep -c "^rs"
done

# VCF statitics

# Indexing the file
for vcf in *_raw_variants_ann.vcf; do
  bgzip -c "$vcf" > "$vcf.gz"
  tabix -p vcf "$vcf.gz"
done

# Calculating the stats
for vcfgz in *_raw_variants_ann.vcf.gz; do
  out=${vcfgz%.vcf.gz}.stats.txt
  JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64 \ # It gave us an error that it needs Java 11 or later, while the Conda environment's Java is 8

    rtg vcfstats "$vcfgz" > "$out"
done

# -------------------------
# Splitting SNPs and INDELs
# -------------------------

for vcf in *_raw_variants_ann.vcf; do
  base=${vcf%.vcf} 

  # SNPs
  gatk --java-options "-Xmx30G" SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V "$vcf" \
    --select-type-to-include SNP \
    -O "${base}_SNP.vcf"

  # INDELs
  gatk --java-options "-Xmx30G" SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V "$vcf" \
    --select-type-to-include INDEL \
    -O "${base}_INDEL.vcf"
done

# -------------------------
# Assessing the different filters in both known and novel variants
# -------------------------

# Working in a separate directory to be organized
cd "$PROJ_DIR/filters" || exit 1

for vcf in ../*_raw_variants_ann_SNP.vcf ../*_raw_variants_ann_INDEL.vcf; do
    base=$(basename "$vcf" .vcf)

    if [[ "$vcf" == *SNP.vcf ]]; then
        var="SNP"
    else
        var="INDEL"
    fi

    for filter in QD MQ MQRankSum FS SOR ReadPosRankSum AN DP InbreedingCoeff; do
        
        out="${base}.${filter}"
        known="known.${base}.${filter}"
        novel="novel.${base}.${filter}"

        awk -v k="${filter}=" '!/^#/{
            n = split($8, a, ";")
            for (i=1; i<=n; i++) {
                if (a[i] ~ "^"k) {
                    sub(k, $3" ", a[i])
                    print a[i]
                }
            }
        }' "$vcf" > "$out"

        grep -v '^\.' "$out" > "$known"
        grep    '^\.' "$out" > "$novel"

    done
done

# -------------------------
# Generating figures
# -------------------------

# Downloading the plotting R script
wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R

for f in HG002_22*_raw_variants_ann_*.*; do
    case "$f" in
        known.*|novel.*) continue ;;
    esac

    label=$(basename "$f") 
    echo "Running Rscript on: $label"
    Rscript densityCurves.R "$label"
done

# -------------------------
# Calculating the DP thresholds
# -------------------------

cat known.HG002_22_*_raw_variants_ann_SNP.DP \
    known.HG002_22_*_raw_variants_ann_INDEL.DP \
  | awk '{
        sum   += $2
        sumsq += ($2)^2
     }
     END {
        mean = sum/NR
        sd   = sqrt((sumsq - sum*sum/NR)/NR)
        print mean, sd, mean + 5*sd
     }'

# 108.212 168.389 950.159

# -------------------------
# SNP Variant filtration
# -------------------------

cd "$PROJ_DIR" || exit 1

for vcf in *_raw_variants_ann_SNP.vcf; do
  base=${vcf%.vcf}

  case "$base" in
    HG002_22_2x_raw_variants_ann_SNP)   SNP_DP_CUTOFF=12  ;;
    HG002_22_10x_raw_variants_ann_SNP)  SNP_DP_CUTOFF=40  ;;
    HG002_22_40x_raw_variants_ann_SNP)  SNP_DP_CUTOFF=110 ;;
    HG002_22_80x_raw_variants_ann_SNP)  SNP_DP_CUTOFF=175 ;;
    HG002_22_raw_variants_ann_SNP)      SNP_DP_CUTOFF=340 ;;
    *) echo "Warning: no DP cutoff defined for $base, skipping"; continue ;;
  esac

  gatk --java-options "-Xmx30G" VariantFiltration \
    -R Homo_sapiens_assembly38.fasta \
    -V "$vcf" \
    --filter-name "snpQD"             --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
    --filter-name "snpMQ"             --filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
    --filter-name "snpMQRankSum"      --filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
    --filter-name "snpFS"             --filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
    --filter-name "snpSOR"            --filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
    --filter-name "snpReadPosRankSum" --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
    --filter-name "snpDP"             --filter-expression "vc.hasAttribute('DP') && DP > ${SNP_DP_CUTOFF}" \
    -O "${base}_clean.vcf"

  echo "Filtered SNP records in ${base}_clean.vcf:"
  grep -v "^#" "${base}_clean.vcf" | awk '$7!="PASS"' | head
done

# -------------------------
# INDEL Variant filtration
# -------------------------

for vcf in *_raw_variants_ann_INDEL.vcf; do
  base=${vcf%.vcf}

  case "$base" in
    HG002_22_2x_raw_variants_ann_INDEL)  INDEL_DP_CUTOFF=12  ;;
    HG002_22_10x_raw_variants_ann_INDEL) INDEL_DP_CUTOFF=40  ;;
    HG002_22_40x_raw_variants_ann_INDEL) INDEL_DP_CUTOFF=110 ;;
    HG002_22_80x_raw_variants_ann_INDEL) INDEL_DP_CUTOFF=175 ;;
    HG002_22_raw_variants_ann_INDEL)     INDEL_DP_CUTOFF=340 ;;
    *) echo "Warning: no DP cutoff defined for $base, skipping"; continue ;;
  esac

  gatk --java-options "-Xmx30G" VariantFiltration \
    -R Homo_sapiens_assembly38.fasta \
    -V "$vcf" \
    --filter-name "indelQD"              --filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
    --filter-name "indelFS"              --filter-expression "vc.hasAttribute('FS') && FS > 200.0" \
    --filter-name "indelSOR"             --filter-expression "vc.hasAttribute('SOR') && SOR > 10.0" \
    --filter-name "indelReadPosRankSum"  --filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
    --filter-name "indelInbreedingCoeff" --filter-expression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
    --filter-name "indelDP"              --filter-expression "vc.hasAttribute('DP') && DP > ${INDEL_DP_CUTOFF}" \
    -O "${base}_clean.vcf"

  echo "Filtered INDEL records in ${base}_clean.vcf:"
  grep -v "^#" "${base}_clean.vcf" | awk '$7!="PASS"' | head
done

# checking the filtered records with a summary report

for vcf in *_clean.vcf; do
  awk -v f="$vcf" '
    BEGIN{t=0;p=0;fl=0;pk=0;pn=0}
    /^#/ {next}
    {t++; if($7=="PASS"){p++; if($3==".") pn++; else pk++} else fl++}
    END{printf "%s\ttotal=%d\tPASS=%d\tfiltered=%d\tPASS_known=%d\tPASS_novel=%d\n", f,t,p,fl,pk,pn}
  ' "$vcf"
done

# Report of how many variants got each filter tag

# SNPs

out="SNP_filter_report.tsv"
: > "$out"

for vcf in *_SNP_clean.vcf; do
  echo "=== $vcf ===" >> "$out"
  echo -e "tag\tfiltered_total\tfiltered_known\tfiltered_novel\t%known_filtered\t%novel_filtered\tknown_total\tnovel_total" >> "$out"

  read KNOWN_TOTAL NOVEL_TOTAL < <(
    awk 'BEGIN{k=0;n=0} /^#/ {next} {if($3==".") n++; else k++} END{print k, n}' "$vcf"
  )

  for tag in snpQD snpMQ snpMQRankSum snpFS snpSOR snpReadPosRankSum snpDP; do
    awk -v t="$tag" -v kt="$KNOWN_TOTAL" -v nt="$NOVEL_TOTAL" '
      BEGIN{tot=0;k=0;n=0}
      /^#/ {next}
      $7!="PASS" {
        m=split($7,a,";")
        hit=0
        for(i=1;i<=m;i++) if(a[i]==t){hit=1; break}
        if(hit){
          tot++
          if($3==".") n++; else k++
        }
      }
      END{
        pk = (kt>0 ? 100*k/kt : 0)
        pn = (nt>0 ? 100*n/nt : 0)
        printf "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%d\n", t, tot, k, n, pk, pn, kt, nt
      }
    ' "$vcf" >> "$out"
  done

  echo "" >> "$out"
done

# INDELs

out="INDEL_filter_report.tsv"
: > "$out"

for vcf in *_INDEL_clean.vcf; do
  echo "=== $vcf ===" >> "$out"
  echo -e "tag\tfiltered_total\tfiltered_known\tfiltered_novel\t%known_filtered\t%novel_filtered\tknown_total\tnovel_total" >> "$out"

  read KNOWN_TOTAL NOVEL_TOTAL < <(
    awk 'BEGIN{k=0;n=0} /^#/ {next} {if($3==".") n++; else k++} END{print k, n}' "$vcf"
  )

  for tag in indelQD indelFS indelSOR indelReadPosRankSum indelInbreedingCoeff indelDP; do
    awk -v t="$tag" -v kt="$KNOWN_TOTAL" -v nt="$NOVEL_TOTAL" '
      BEGIN{tot=0;k=0;n=0}
      /^#/ {next}
      $7!="PASS" {
        m=split($7,a,";")
        hit=0
        for(i=1;i<=m;i++) if(a[i]==t){hit=1; break}
        if(hit){
          tot++
          if($3==".") n++; else k++
        }
      }
      END{
        pk = (kt>0 ? 100*k/kt : 0)
        pn = (nt>0 ? 100*n/nt : 0)
        printf "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%d\n", t, tot, k, n, pk, pn, kt, nt
      }
    ' "$vcf" >> "$out"
  done

  echo "" >> "$out"
done

# -------------------------
# Starting benchmarking
# -------------------------

cd "$PROJ_DIR" || exit 1

# -------------------------
# Downloading GIAB HG002 (GRCh38) truth VCF and confident regions BED, and preparing them for the analysis
# -------------------------

TRUTH_VCF_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

wget -O benchmarking/truth/HG002_truth_GRCh38_v4.2.1.vcf.gz "$TRUTH_VCF_URL"
wget -O benchmarking/bed/HG002_confident_GRCh38_v4.2.1.bed "$TRUTH_BED_URL"

# Indexing the full truth VCF
tabix -p vcf benchmarking/truth/HG002_truth_GRCh38_v4.2.1.vcf.gz

# Extracting confident regions restricted to chr22

awk '$1=="chr22"' benchmarking/bed/HG002_confident_GRCh38_v4.2.1.bed \
  > benchmarking/bed/HG002_confident_chr22.bed

# Sorting chr22 confident regions
sort -k1,1 -k2,2n benchmarking/bed/HG002_confident_chr22.bed \
  -o benchmarking/bed/HG002_confident_chr22.bed

# Compressing and indexing chr22 confident regions BED
bgzip -f benchmarking/bed/HG002_confident_chr22.bed
tabix -p bed benchmarking/bed/HG002_confident_chr22.bed.gz

# Extracting chr22 variants from the truth VCF
bcftools view -Oz -r chr22 \
  -o benchmarking/truth/HG002_truth_chr22.vcf.gz \
  benchmarking/truth/HG002_truth_GRCh38_v4.2.1.vcf.gz

# Indexing chr22 truth VCF
tabix -p vcf benchmarking/truth/HG002_truth_chr22.vcf.gz

# Splitting chr22 truth VCF into SNPs and INDELs

# SNPS
bcftools view -v snps -Oz \
  -o benchmarking/truth/HG002_truth_chr22_SNP.vcf.gz \
  benchmarking/truth/HG002_truth_chr22.vcf.gz

# INDELs
bcftools view -v indels -Oz \
  -o benchmarking/truth/HG002_truth_chr22_INDEL.vcf.gz \
  benchmarking/truth/HG002_truth_chr22.vcf.gz

# Indexing SNPs and INDELs truth VCFs
tabix -p vcf benchmarking/truth/HG002_truth_chr22_SNP.vcf.gz
tabix -p vcf benchmarking/truth/HG002_truth_chr22_INDEL.vcf.gz

# Extracting chr22 regions from the exome capture BED

awk '$1=="chr22"' S04380110_Regions.bed \
  | sort -k1,1 -k2,2n \
  > benchmarking/bed/exome_chr22.bed

# # Restricting confident regions to exome capture targets (intersecting GIAB confident chr22 with exome chr22)

bedtools intersect \
  -a benchmarking/bed/HG002_confident_chr22.bed.gz \
  -b benchmarking/bed/exome_chr22.bed \
  > benchmarking/bed/HG002_confident_exome_chr22.bed

# Sorting intersected confident exome BED

sort -k1,1 -k2,2n benchmarking/bed/HG002_confident_exome_chr22.bed \
  -o benchmarking/bed/HG002_confident_exome_chr22.bed

# Compressing and indexing confident exome chr22 BED

bgzip -f benchmarking/bed/HG002_confident_exome_chr22.bed
tabix -p bed benchmarking/bed/HG002_confident_exome_chr22.bed.gz

# -------------------------
# hap.py
# -------------------------

# Creating a dedicated Conda environment for hap.py
# hap.py installation requires Python 2.7, which is incompatible with the ngs1 environment

conda create -n happyenv python=2.7
conda activate happyenv

# Creating RTG SDF reference for hap.py vcfeval

rtg format \
  -o Homo_sapiens_assembly38.sdf \
  Homo_sapiens_assembly38.fasta

# Defining key paths for hap.py and vcfeval

export REF="$PWD/Homo_sapiens_assembly38.fasta"
export SDF="$PWD/Homo_sapiens_assembly38.sdf"
export TRUTH_SNP="$PWD/benchmarking/truth/HG002_truth_chr22_SNP.vcf.gz"
export TRUTH_INDEL="$PWD/benchmarking/truth/HG002_truth_chr22_INDEL.vcf.gz"
export BED="$PWD/benchmarking/bed/HG002_confident_exome_chr22.bed.gz"

# Running hap.py per coverage and per variant type

cd "$PROJ_DIR/benchmarking/results" || exit 1

for cov in 2x 10x 40x 80x; do
  echo "=== $(date) | $cov SNP ==="
  hap.py "$TRUTH_SNP" "../../HG002_22_${cov}_raw_variants_ann_SNP_clean.vcf.gz" \
    -f "$BED" \
    -r "$REF" \
    --engine vcfeval \
    --engine-vcfeval-template "$SDF" \
    --threads 14 \
    --pass-only \
    -o "HG002_chr22_${cov}_SNP"

  echo "=== $(date) | $cov INDEL ==="
  hap.py "$TRUTH_INDEL" "../../HG002_22_${cov}_raw_variants_ann_INDEL_clean.vcf.gz" \
    -f "$BED" \
    -r "$REF" \
    --engine vcfeval \
    --engine-vcfeval-template "$SDF" \
    --threads 14 \
    --pass-only \
    -o "HG002_chr22_${cov}_INDEL"
done
# -------------------------
# Summarizing hap.py results (PASS only from hap.py outputs)
# -------------------------

EVAL_BP=$(
  zcat -f "$BED" \
  | awk 'BEGIN{sum=0} {sum+=($3-$2)} END{print sum}'
)

{
  echo -e "cov\ttype\tfilter\ttruth_total\ttruth_tp\ttruth_fn\tquery_total\tquery_fp\tquery_unk\trecall\tprecision\tf1\tgt_concordance\tfp_per_mb"

  for cov in 2x 10x 40x 80x; do
    for type in SNP INDEL; do

      sum="HG002_chr22_${cov}_${type}.summary.csv"
      ext="HG002_chr22_${cov}_${type}.extended.csv"

      GT_CONC=$(
        awk -F',' -v t="$type" '
          NR==1{
            for(i=1;i<=NF;i++){
              if($i ~ /(GT.*Concord|Genotype.*Concord|Concordance.*GT)/){c=i}
            }
            next
          }
          ($1==t && $2=="PASS"){print (c? $c : "NA"); exit}
        ' "$ext"
      )

      awk -F',' -v cov="$cov" -v type="$type" -v gt="$GT_CONC" -v bp="$EVAL_BP" '
        NR==1 { for (i=1; i<=NF; i++) h[$i]=i; next }
        ($1==type && $2=="PASS") {
          fp = $(h["QUERY.FP"]) + 0
          fpm = (bp>0 ? fp/(bp/1000000.0) : "NA")

          printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            cov, $1, $2,
            $(h["TRUTH.TOTAL"]), $(h["TRUTH.TP"]), $(h["TRUTH.FN"]),
            $(h["QUERY.TOTAL"]), $(h["QUERY.FP"]), $(h["QUERY.UNK"]),
            $(h["METRIC.Recall"]), $(h["METRIC.Precision"]), $(h["METRIC.F1_Score"]),
            gt, fpm
        }
      ' "$sum"

    done
  done
} > benchmarking_summary_chr22.tsv

# -------------------------
# Plotting benchmarking metrics
# -------------------------

# Switching back to the main environment for plotting
conda activate ngs1

# Directory setting 
cd "$PROJ_DIR" || exit 1

Rscript - <<'RS'
df <- read.table("benchmarking/results/benchmarking_summary_chr22.tsv",
                 header=TRUE, sep="\t", stringsAsFactors=FALSE)

df$cov_num <- as.numeric(sub("x","", df$cov))

plot_metric <- function(metric, ylab, outpng, ylim=NULL) {
  png(outpng, width=1200, height=800)

  if (is.null(ylim)) {
    ylim <- range(df[[metric]], na.rm=TRUE)
    if (all(is.finite(ylim)) && ylim[1] == ylim[2]) ylim <- c(0, ylim[2])
  }

  plot(NA, xlim=range(df$cov_num), ylim=ylim,
       xlab="Coverage (x)", ylab=ylab,
       main=paste0(ylab, " vs Coverage (chr22, confident & exome)"))

  for (t in c("SNP","INDEL")) {
    d <- df[df$type==t, ]
    d <- d[order(d$cov_num), ]

    lines(d$cov_num, d[[metric]], type="b")
    text(d$cov_num, d[[metric]], labels=t, pos=3, cex=0.8)
  }

  dev.off()
}

plot_metric("recall",          "Recall",                 "benchmarking/results/plot_recall.png",    ylim=c(0,1))
plot_metric("precision",       "Precision",              "benchmarking/results/plot_precision.png", ylim=c(0,1))
plot_metric("f1",              "F1 score",               "benchmarking/results/plot_f1.png",        ylim=c(0,1))
plot_metric("gt_concordance",  "Genotype concordance",   "benchmarking/results/plot_gt_conc.png",   ylim=c(0,1))
plot_metric("fp_per_mb",       "False positives per Mb", "benchmarking/results/plot_fp_per_mb.png")
RS
