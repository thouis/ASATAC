# TODO:
# adapter clipping - fastq-mcf
# read group assignment - from metadata?
# CleanSam - needed after bwa mem?


# Reference data fetched via "wget -r ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/"
# Then unzip ucsc.hg19 and build its index with bwa index.
genome = "reference/ucsc.hg19.fasta"
known_vcf = "reference/dbsnp_138.hg19.vcf.gz"

# 00-common from 1K genomes, filtered to SNPs with single alternates
count_vcf = "other_vcf/00-common-unique-biallelic.vcf"

# paths to java-based tools
PICARD = "java -jar ~/bin/picard.jar"
GATK = "java -jar ~/bin/GenomeAnalysisTK.jar"

rule all:
    input:
        "allele_counts/SRR2378583_subset.counts.csv"

rule subset_fastq:
    input:
        "fastq/{sample}_{index}.fastq.gz"
    output:
        "fastq/{sample}_subset_{index}.fastq.gz"
    shell:
        "zcat {input} | sed -n 1,4000000p | gzip -c > {output}"

rule bwa_map:
    input:
        "fastq/{sample}_1.fastq.gz",
        "fastq/{sample}_2.fastq.gz"
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg = "@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tPU:1\\tLB:lib1"
    shell:
        # rg should be in a config file (one per SRA/fastq, probably)
        # adapter trimming probably should be, as well.
        "bwa mem -M -t 16 -R '{params.rg}' {genome} {input} | samtools view -Sbh - > {output}"

rule samtools_index:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "mapped_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bam_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("mapped_reads/{sample}.sorted.bam")
    shell:
        "samtools sort -o {input} _ > {output}"

rule mark_duplicates:
    input:
        "mapped_reads/{sample}.sorted.bam"
    output:
        bam="mapped_reads/{sample}.sorted.dedup.bam"
        metrics="metrics/markduplicates_{sample}.txt"
    shell:
        ("{PICARD} MarkDuplicates"
         " CREATE_INDEX=true TMP_DIR=/tmp"
         " INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics}")

rule BQSR:
    input:
        bam="mapped_reads/{sample}.dedup.bam",
        bai="mapped_reads/{sample}.dedup.bam.bai"
    output:
        table="metrics/calibration_{sample}.table",
        bam="mapped_reads/{sample}.dedup.recal.bam"
    params:
        br = ("--default_platform illumina -cov QualityScoreCovariate -cov ContextCovariate"
              " --downsampling_type NONE -nct 4"),
        pr = "-nct 4"
    shell:
        """
        {GATK} -T BaseRecalibrator {params.br} -R {genome}  -knownSites {known_vcf} -I {input.bam} -o {output.table}
        {GATK} -T PrintReads {params.pr} -I {input.bam} -BQSR {output.table} -R {genome} -o {output.bam}
        """
        # do indel realignment?

rule counts:
    input:
        bam="mapped_reads/{sample}.sorted.dedup.recal.bam",
        bai="mapped_reads/{sample}.sorted.dedup.recal.bam.bai",
    output:
        "allele_counts/{sample}.counts.csv"
    shell:
        "{GATK} -T ASEReadCounter -R {genome} -o {output} -I {input.bam} -sites {count_vcf}"

# useful sources for the rules above.
# https://nextflow-io.github.io/hack17-varcall/
# https://gatkforums.broadinstitute.org/gatk/discussion/9622/allele-specific-annotation-and-filtering
# http://www.cureffi.org/2013/08/26/allele-specific-rna-seq-pipeline-using-gsnap-and-gatk/
# http://barcwiki.wi.mit.edu/wiki/SOPs/variant_calling_GATK
