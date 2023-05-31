#!/bin/bash

sample=$1
ref_length_bed=$2
ref_bed=$3

# shellcheck disable=code

# This code is used to extract the statistics of clean data from json file
# variables: raw_bases, raw_reads, trimmed_bases, trimmed_reads, Q20_bases, Q30_bases
raw_bases=`grep "total_bases" 01.clean_data/${sample}_fastp.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`
raw_reads=`grep "total_reads" 01.clean_data/${sample}_fastp.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`
trimmed_bases=`grep "total_bases" 01.clean_data/${sample}_fastp.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n2 | tail -n1`
trimmed_reads=`grep "total_reads" 01.clean_data/${sample}_fastp.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n2 | tail -n1`
Q20_bases=`grep "q20_bases" 01.clean_data/${sample}_fastp.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`
Q30_bases=`grep "q30_bases" 01.clean_data/${sample}_fastp.json | awk 'BEGIN{FS=":|: |,"}{print $2}' | head -n1`

# Get the number of reads that mapped to the genome
mapped_reads_includes=`samtools flagstat 02.bismark_map/${sample}_clean_1_bismark_bt2_pe.bam | grep "mapped" | awk 'BEGIN{FS=" "}{print $1}' | head -n1`

# Get the number of reads that mapped to the genome as a secondary alignment
mapped_reads_secondary=`samtools flagstat 02.bismark_map/${sample}_clean_1_bismark_bt2_pe.bam | grep "secondary" | awk 'BEGIN{FS=" "}{print $1}' | head -n1`

# Get the number of reads that mapped to the genome as a supplementary alignment
mapped_reads_supplementary=`samtools flagstat 02.bismark_map/${sample}_clean_1_bismark_bt2_pe.bam | grep "supplementary" | awk 'BEGIN{FS=" "}{print $1}' | head -n1`

# Calculate the number of reads that mapped to the genome
mapped_reads=$((${mapped_reads_includes}-${mapped_reads_secondary}-${mapped_reads_supplementary}))

# Calculate the percentage of reads that mapped to the genome
mapped_ratio=$(printf "%.2f" `echo "scale=2;100*${mapped_reads}/${trimmed_reads}"|bc`)

# This code counts the number of reads that were mapped in the deduplication report, and multiplies it by two since we have paired-end reads. This number is used to calculate the percentage of reads that were removed during deduplication in the next code chunk.
dup_read_pairs=$(cat 02.bismark_map/${sample}.deduplication_report.txt | cut -f 2 | head -n 3 | tail -n 1 | cut -d " " -f 1)
dup_reads=$(($dup_read_pairs*2))

# Calculate the number of deduplicated reads and the deduplication ratio
dedup_reads=$((${mapped_reads}-${dup_reads}))
dedup_ratio=$(printf "%.2f" `echo "scale=2;100-100*${dedup_reads}/${mapped_reads}"|bc`)

# 修改coverage、ratio、depth等函数
# coverage=`bedtools merge -i 02.map/${sample}_mapQ30_sort_rmdup.bed | bedtools intersect -a - -b ${panel_bed} | sort -V -k1,1 -k2,2 -i - | bedtools jaccard -a - -b ${panel_bed} | awk 'BEGIN{FS="\t"}{print $3}' | tail -n1`
# 总的覆盖区域长度
coverage=`bedtools merge -i 02.bismark_map/${sample}_sorted_rmdup.bed | bedtools sort -i - | bedtools jaccard -a ${ref_bed} -b - | awk 'BEGIN{FS="\t"}{print $3}' | tail -n1`
# 在Genome上的深度
depth=`bedtools intersect -a ${ref_bed} -b 02.bismark_map/${sample}_sorted_rmdup.bed | bedtools genomecov -i - -g ${ref_length_bed} | awk '{FS=OFS="\t"}{ if ($1=="genome"){sum_i+=$2*$3;total_len=$4}}END{print sum_i/total_len}'`

# 甲基化覆盖率
lambda_reads=`grep -c "lambda" 02.bismark_map/${sample}_sorted_rmdup.bed`
conversion=`zcat 03.extracted_met_sites/${sample}.deduplicated.CX_report.txt.gz | grep "lambda" | awk 'BEGIN{FS="\t"}{methy+=$4;total+=$4+$5}END{print (1-methy/total)}'`
genome_CG=`zcat 03.extracted_met_sites/${sample}.deduplicated.CX_report.txt.gz | grep "chr" | awk 'BEGIN{FS="\t"}{ if ($6=="CG") {methy+=$4;total+=$4+$5}}END{print methy/total}'`
genome_CHH=`zcat 03.extracted_met_sites/${sample}.deduplicated.CX_report.txt.gz | grep "chr" | awk 'BEGIN{FS="\t"}{ if ($6=="CHH") {methy+=$4;total+=$4+$5}}END{print methy/total}'`
genome_CHG=`zcat 03.extracted_met_sites/${sample}.deduplicated.CX_report.txt.gz | grep "chr" | awk 'BEGIN{FS="\t"}{ if ($6=="CHG") {methy+=$4;total+=$4+$5}}END{print methy/total}'`

echo -e sample'\t'raw_bases'\t'raw_reads'\t'Q20_bases'\t'Q30_bases'\t'trimmed_bases'\t'trimmed_reads'\t'mapped_reads'\t'mapped_ratio'\t'dedup_reads'\t'dedup_ratio'\t'coverage'\t'depth'\t'lambda_reads'\t'conversion'\t'genome_CG'\t'genome_CHH'\t'genome_CHG > 04.stat_info/${sample}_met_stats.tsv
echo -e ${sample}'\t'${raw_bases}'\t'${raw_reads}'\t'${Q20_bases}'\t'${Q30_bases}'\t'${trimmed_bases}'\t'${trimmed_reads}'\t'${mapped_reads}'\t'${mapped_ratio}'\t'${dedup_reads}'\t'${dedup_ratio}'\t'${coverage}'\t'${depth}'\t'${lambda_reads}'\t'${conversion}'\t'${genome_CG}'\t'${genome_CHH}'\t'${genome_CHG} >> 04.stat_info/${sample}_met_stats.tsv

echo -e "["$(date)"]\t$sample statistics done!"