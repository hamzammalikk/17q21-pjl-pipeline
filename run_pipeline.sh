#!/usr/bin/env bash
# 17q21 asthma locus pipeline: 1000 Genomes PJL, GRCh38, window chr17:39700000-40000000
# Reproduces the steps that were run interactively. Run inside Ubuntu (WSL2):  bash run_pipeline.sh
#
# Needs: the conda env "pjl17q21" (samtools, bcftools, bwa-mem2, fastp, plink, plink2, python with pandas/numpy/matplotlib)
# Needs these helper scripts in ~/pjl_17q21/scripts/ (see the handoff document): ld_heatmap.py, vep_rest.py, eqtl_lookup.py, merge_final.py
# Most download/heavy steps are skipped if their output already exists, so re-running is safe.

set -euo pipefail

mkdir -p ~/pjl_17q21
cat > ~/pjl_17q21/env.sh << 'EOF'
export PROJ=$HOME/pjl_17q21
export REGION=chr17:39700000-40000000
export BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
conda activate pjl17q21
cd $PROJ
EOF
# shellcheck disable=SC1090
source "$(conda info --base)/etc/profile.d/conda.sh"
source ~/pjl_17q21/env.sh
mkdir -p ref data/reads results/bam results/qc logs scripts

S3=https://1000genomes.s3.amazonaws.com/1000G_2504_high_coverage/data
LEAD_ID=17:39913696:C:T          # rs7216389 (GRCh38 chr17:39913696 C/T), the lead SNP used for LD
step() { echo; echo "=== [$(date +%H:%M:%S)] $*"; }

# ---------------------------------------------------------------- Module 1
step "1. PJL sample list (96 unrelated) and regional VCF"
[ -s 20130606_g1k_3202_samples_ped_population.txt ] || curl -sO "$BASE/20130606_g1k_3202_samples_ped_population.txt"
[ -s 1000G_2504_high_coverage.sequence.index ]      || curl -sO "$BASE/1000G_2504_high_coverage.sequence.index"
awk '$6=="PJL"{print $2}' 20130606_g1k_3202_samples_ped_population.txt > data/pjl.all.txt      # 146 incl. relatives
grep -o -w -F -f data/pjl.all.txt 1000G_2504_high_coverage.sequence.index | sort -u > data/pjl.samples.txt
echo "unrelated PJL samples: $(wc -l < data/pjl.samples.txt) (expect 96)"

VCF=$BASE/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
if [ ! -s data/pjl_17q21.vcf.gz.tbi ]; then
  bcftools view -r "$REGION" -S data/pjl.samples.txt "$VCF" -Oz -o data/pjl_17q21.vcf.gz
  bcftools index -t data/pjl_17q21.vcf.gz
fi
echo "regional VCF: $(bcftools query -l data/pjl_17q21.vcf.gz | wc -l) samples, $(bcftools index -n data/pjl_17q21.vcf.gz) records"

# ---------------------------------------------------------------- Module 2
step "2a. chr17 reference (Ensembl, header renamed to chr17) + bwa-mem2 index"
if [ ! -s ref/chr17.fa.bwt.2bit.64 ]; then
  curl -C - --retry 10 --retry-all-errors -o ref/chr17.raw.fa.gz \
    https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz
  gzip -t ref/chr17.raw.fa.gz
  gunzip -c ref/chr17.raw.fa.gz | sed '1s/.*/>chr17/' > ref/chr17.fa
  samtools faidx ref/chr17.fa
  bwa-mem2 index ref/chr17.fa
fi
cut -f1,2 ref/chr17.fa.fai        # expect: chr17  83257441

step "2b. reads -> FASTQ -> fastp -> bwa-mem2 -> sorted, duplicate-marked BAM (5 samples)"
head -5 data/pjl.samples.txt > data/aln.samples.txt
for S in $(cat data/aln.samples.txt); do
  if [ -s "results/bam/$S.bam.bai" ]; then echo "$S: BAM exists, skipping"; continue; fi
  ERR=$(grep -w "$S" 1000G_2504_high_coverage.sequence.index | grep -o 'ERR[0-9]\{7\}' | head -1)
  [ -n "$ERR" ] || { echo "no run id found for $S"; exit 1; }
  echo "$S $ERR"
  samtools view -T ref/chr17.fa -b -F 0x900 -o "data/reads/$S.region.bam" "$S3/$ERR/$S.final.cram" "$REGION"
  samtools sort -n -@2 -o "data/reads/$S.nsort.bam" "data/reads/$S.region.bam"
  samtools fastq -@2 -n -1 "data/reads/${S}_R1.fq.gz" -2 "data/reads/${S}_R2.fq.gz" -0 /dev/null -s /dev/null "data/reads/$S.nsort.bam"
  fastp -i "data/reads/${S}_R1.fq.gz" -I "data/reads/${S}_R2.fq.gz" \
        -o "data/reads/${S}_R1.clean.fq.gz" -O "data/reads/${S}_R2.clean.fq.gz" \
        --detect_adapter_for_pe -q 20 -l 50 -w 4 -j "results/qc/$S.fastp.json" -h "results/qc/$S.fastp.html"
  bwa-mem2 mem -t 4 -R "@RG\tID:$S\tSM:$S\tPL:ILLUMINA\tLB:$S" ref/chr17.fa \
      "data/reads/${S}_R1.clean.fq.gz" "data/reads/${S}_R2.clean.fq.gz" \
   | samtools fixmate -m - - \
   | samtools sort -@4 -T "results/bam/tmp_$S" - \
   | samtools markdup - "results/bam/$S.bam"
  samtools index "results/bam/$S.bam"
done
for S in $(cat data/aln.samples.txt); do
  echo "== $S"; samtools flagstat "results/bam/$S.bam" | grep -E "in total|mapped \(|properly paired"
  samtools coverage -r "$REGION" "results/bam/$S.bam" | cut -f1-9
done

# ---------------------------------------------------------------- Module 3
step "3. variant calling, filtering, benchmark vs 1000G"
ls results/bam/*.bam > results/bamlist.txt
bcftools mpileup -f ref/chr17.fa -r "$REGION" -q 20 -Q 20 -a FORMAT/AD,FORMAT/DP -b results/bamlist.txt -Ou \
  | bcftools call -mv -Oz -o results/calls.raw.vcf.gz
bcftools index -t results/calls.raw.vcf.gz
bcftools norm -f ref/chr17.fa -m -any results/calls.raw.vcf.gz -Ou \
  | bcftools filter -e 'QUAL<30' -Ou \
  | bcftools filter -S . -e 'FMT/DP<8' -Oz -o results/calls.filt.vcf.gz
bcftools index -t results/calls.filt.vcf.gz
bcftools stats results/calls.filt.vcf.gz | grep -E "number of (SNPs|indels):|^TSTV"

SAMPLES=$(paste -sd, data/aln.samples.txt)
bcftools view -s "$SAMPLES" -c1 -v snps -m2 -M2 data/pjl_17q21.vcf.gz -Oz -o results/truth.snps.vcf.gz
bcftools view -v snps -m2 -M2 results/calls.filt.vcf.gz -Oz -o results/calls.snps.vcf.gz
bcftools index -t results/truth.snps.vcf.gz
bcftools index -t results/calls.snps.vcf.gz
rm -rf results/isec
bcftools isec -p results/isec -Oz results/calls.snps.vcf.gz results/truth.snps.vcf.gz
for f in results/isec/000{0,1,2,3}.vcf.gz; do echo "$f $(bcftools view -H "$f" | wc -l)"; done
echo "(0000 = only ours, 0001 = only 1000G, 0002/0003 = shared)"

# ---------------------------------------------------------------- Module 4
step "4. PLINK conversion, QC, LD with the lead SNP, haplotype blocks, LD heatmap"
plink2 --vcf data/pjl_17q21.vcf.gz --double-id --snps-only just-acgt --max-alleles 2 \
       --set-all-var-ids '@:#:$r:$a' --make-bed --out results/pjl17
plink2 --bfile results/pjl17 --maf 0.01 --geno 0.05 --hwe 1e-6 --make-bed --out results/pjl17.qc
plink2 --bfile results/pjl17.qc --freq --out results/pjl17.qc
echo "SNPs after QC: $(wc -l < results/pjl17.qc.bim) (expect 850)"
plink --bfile results/pjl17.qc --r2 --ld-snp "$LEAD_ID" --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0 --out results/ld_lead
plink --bfile results/pjl17.qc --blocks no-pheno-req --blocks-max-kb 500 --out results/blocks
plink2 --bfile results/pjl17.qc --maf 0.10 --make-bed --out results/ld_sub
plink --bfile results/ld_sub --r2 square --out results/ld_matrix
python3 scripts/ld_heatmap.py

# ---------------------------------------------------------------- Module 5
step "5. candidates (r2 >= 0.6 with lead SNP), VEP (Ensembl REST), GTEx eQTLs, final table"
awk 'NR>1 && $7>=0.6 {print "chr17\t"$5} NR==2 {print "chr17\t"$2}' results/ld_lead.ld | sort -k2,2n -u > results/cand.pos.tsv
bcftools view -T results/cand.pos.tsv -G -v snps data/pjl_17q21.vcf.gz -Oz -o results/cand.sites.vcf.gz
bcftools index -t results/cand.sites.vcf.gz
python3 scripts/vep_rest.py
python3 scripts/eqtl_lookup.py
python3 scripts/merge_final.py

step "DONE. Main outputs: results/final_table.tsv, results/ld_heatmap.png, results/ld_lead.ld, results/blocks.blocks.det, results/vep.tsv, results/eqtl_hits.tsv"
