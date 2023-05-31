import os 

file_list = os.listdir("00.raw_data")
name_list = set([x.split('.')[0][:-2] for x in file_list])
gatk_ref_path = "/datb1/wuyuxuan/resources/GATK4_hg38" # includes GATK4_hg38 data
bismark_ref_path = "/datb1/wuyuxuan/resources/bismark_hg38_with_lambda" # includes bismark data

"""
onsuccess:
    shell("python /home/wuyuxuan/py_scripts/notify.py")

onerror:
    shell("python /home/wuyuxuan/py_scripts/notify.py Error! an error occured")
"""

rule all:
    input:
        cgsites_met_bz=expand("03.extracted_met_sites/{sample}_CG_sites.bed.gz", sample=name_list),
        stats=expand("04.stat_info/{sample}_met_stats.tsv", sample=name_list),
        atac_narrowpeak=expand("05.genrich/{sample}.narrowPeak", sample=name_list),
        

## step1 - filter
rule TrimData:
    input:
        R1="00.raw_data/{sample}_1.fq.gz",
        R2="00.raw_data/{sample}_2.fq.gz"
    output:
        clean_R1=temp("01.clean_data/{sample}_clean_1.fq.gz"),
        clean_R2=temp("01.clean_data/{sample}_clean_2.fq.gz"),
        report_json="01.clean_data/{sample}_fastp.json",
        report_html="01.clean_data/{sample}_fastp.html"
    threads: 4
    conda:
        "ngspipe"
    shell: """
        fastp -D -i {input.R1} -I {input.R2} -o {output.clean_R1} -O {output.clean_R2} \
        -j {output.report_json} -h {output.report_html} 
        echo -e "["$(date)"]\\t{wildcards.sample} QC done!"
    """

rule Multiqc_fastp:
    input:
        report_json=expand("01.clean_data/{sample}_fastp.json", sample=name_list),
        report_html=expand("01.clean_data/{sample}_fastp.html", sample=name_list)
    output:
        multiqc_report="01.clean_data/multiqc_fastp.html"
    threads: 1
    conda: "tgspipe"
    shell: """
        multiqc -f -o 01.clean_data -n multiqc_fastp 01.clean_data
    """

## step2 - bismark map
rule BS_mapping:
    input:
        clean_R1="01.clean_data/{sample}_clean_1.fq.gz",
        clean_R2="01.clean_data/{sample}_clean_2.fq.gz",
    output: 
        bis_mapped_bam="02.bismark_map/{sample}_clean_1_bismark_bt2_pe.bam",
        bis_rmdup_bam="02.bismark_map/{sample}.deduplicated.bam",
        bis_rmdup_report="02.bismark_map/{sample}.deduplication_report.txt",
        bis_mapping_report="02.bismark_map/{sample}_clean_1_bismark_bt2_PE_report.txt",
    threads: 4
    conda:
        "bs_map"
    shadow: "full"
    shell: """
        bismark -p {threads} --parallel 2 --fastq --non_directional --output_dir ./02.bismark_map {bismark_ref_path} -1 {input.clean_R1} -2 {input.clean_R2}
        deduplicate_bismark -p --output_dir ./02.bismark_map -o {wildcards.sample} --bam {output.bis_mapped_bam}
    """

## step3 - extract methylation
rule ExtractMethyl:
    input: 
        bis_rmdup_bam="02.bismark_map/{sample}.deduplicated.bam",
    output:
        cgsites_met_bz="03.extracted_met_sites/{sample}_CG_sites.bed.gz",
        chhsites_met_bz="03.extracted_met_sites/{sample}_CHH_sites.bed.gz",
        chgsites_met_bz="03.extracted_met_sites/{sample}_CHG_sites.bed.gz",
        CX_report="03.extracted_met_sites/{sample}.deduplicated.CX_report.txt.gz",
    threads: 4
    conda:
        "bs_map"
    shadow: "full"
    shell: """
    bismark_methylation_extractor --cytosine_report --CX --gzip --multicore {threads} --o ./03.extracted_met_sites --buffer_size 8G --genome_folder {bismark_ref_path} {input.bis_rmdup_bam}
    zcat {output.CX_report} | grep "chr" | awk 'BEGIN {{FS=OFS="\\t"}}{{ if (($6=="CG") && ($4+$5!=0)) {{print $1,$2,$2,$6,$4/($4+$5),$3}} }}' | bgzip > {output.cgsites_met_bz}
    tabix -f {output.cgsites_met_bz}
    zcat {output.CX_report} | grep "chr" | awk 'BEGIN {{FS=OFS="\\t"}}{{ if (($6=="CHH") && ($4+$5!=0)) {{print $1,$2,$2,$6,$4/($4+$5),$3}} }}' | bgzip > {output.chhsites_met_bz}
    tabix -f {output.chhsites_met_bz}
    zcat {output.CX_report} | grep "chr" | awk 'BEGIN {{FS=OFS="\\t"}}{{ if (($6=="CHG") && ($4+$5!=0)) {{print $1,$2,$2,$6,$4/($4+$5),$3}} }}' | bgzip > {output.chgsites_met_bz}
    tabix -f {output.chgsites_met_bz}
    """
# | grep -E "chr[0-9]{1,2}\s|chr[XY]" capture chr1-22|XY

## step4 - stats info
rule Stats:
    input:
        report_json="01.clean_data/{sample}_fastp.json",
        report_html="01.clean_data/{sample}_fastp.html",
        bis_mapped_bam="02.bismark_map/{sample}_clean_1_bismark_bt2_pe.bam",
        bis_rmdup_bam="02.bismark_map/{sample}.deduplicated.bam",
        bis_rmdup_report="02.bismark_map/{sample}.deduplication_report.txt",
        CX_report="03.extracted_met_sites/{sample}.deduplicated.CX_report.txt.gz",
        chrom_sizes="ref/hg38.chrom.sizes",
        ref_bed="ref/hg38.genome.bed",
        script="stats_met_genome.sh",
    output:
        stats="04.stat_info/{sample}_met_stats.tsv",
        sorted_rmdup_bam="02.bismark_map/{sample}_sorted_rmdup.bam",
        sorted_rmdup_bed="02.bismark_map/{sample}_sorted_rmdup.bed",
    threads: 4
    conda:
        "ngspipe"
    shell: """
        samtools sort -@ {threads} {input.bis_rmdup_bam} -o {output.sorted_rmdup_bam}
        bedtools bamtobed -i {output.sorted_rmdup_bam} > {output.sorted_rmdup_bed}
        sh {input.script} {wildcards.sample} {input.chrom_sizes} {input.ref_bed}
    """

## step5 - genrich peaks calling
rule Genrich:
    input:
        sorted_rmdup_bam="02.bismark_map/{sample}_sorted_rmdup.bam",
    output:
        no_mt_bam=temp("05.genrich/{sample}_no_chrMT.bam"),
        insert_size_report="05.genrich/{sample}_insert_size.pdf",
        no_mt_sort_bam="05.genrich/{sample}_no_chrMT_sort.bam",
        atac_narrowpeak="05.genrich/{sample}.narrowPeak",
        atac_bed="05.genrich/{sample}.bed",
        frip_stats="05.genrich/{sample}_frip.txt",
    threads: 4  
    conda:
        "ngspipe"
    shell: """
        samtools index {input.sorted_rmdup_bam}
        samtools idxstats {input.sorted_rmdup_bam} | cut -f 1 | grep -E "chr[0-9]{{1,2}}|chr[XY]" |\
        xargs samtools view -o {output.no_mt_bam} {input.sorted_rmdup_bam}
        conda activate variant_calling
        samtools index {output.no_mt_bam}
        gatk CollectInsertSizeMetrics -H {output.insert_size_report} -I {output.no_mt_bam} -O 05.genrich/{wildcards.sample}_insert_size.txt
        samtools sort -n -@ 4 -m 200000000 -o {output.no_mt_sort_bam} {output.no_mt_bam}
        conda activate atac_seq
        Genrich -t {output.no_mt_sort_bam} -o {output.atac_narrowpeak} -b {output.atac_bed} -m 30 -q 0.05 -r -j
        reads_total=$(samtools view -c {output.no_mt_sort_bam})
        reads_inregion=$(samtools view -c -L {output.atac_narrowpeak} {output.no_mt_sort_bam})
        printf "{wildcards.sample} frip: %.2f\n" `echo "scale=2;$reads_inregion/$reads_total"|bc` > {output.frip_stats}
    """


