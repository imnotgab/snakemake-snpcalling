SAMPLES = ["TLE66_N", "TLE66_T"]
DB = "/staging/leuven/stg_00079/teaching/hg38_9/chr9.fa"
SNPEFF_JAR = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2026/share/snpeff-5.4.0a-0/snpEff.jar"
SNPEFF_DB = "/staging/leuven/stg_00079/teaching/snpeff_db"

rule all:
    input:
        "100.final/snps.annotated.tsv",
        expand("010.fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("020.bam/{sample}.flagstat", sample=SAMPLES)

rule fastqc:
    input:
        reads="{sample}.fastq"
    output:
        html="010.fastqc/{sample}_fastqc.html",
        zip="010.fastqc/{sample}_fastqc.zip"
    shell:
        """
        test -s {input.reads}
        head -n 1 {input.reads} | grep -q "^@"

        mkdir -p 010.fastqc
        fastqc {input.reads} -o 010.fastqc/

        echo "testing fastqc..."
        test -s {output.html}
        """

rule bwa_mem_and_sort:
    input:
        reads="{sample}.fastq"
    output:
        bam="020.bam/{sample}.bam"
    shell:
        """
        test -s {input.reads}
        head -n 1 {input.reads} | grep -q "^@"

        mkdir -p 020.bam
        bwa mem {DB} {input.reads} | samtools sort -o {output.bam} -
        
        echo "testing pliku BAM..."
        samtools quickcheck {output.bam}
        """

rule samtools_index:
    input:
        bam="020.bam/{sample}.bam"
    output:
        bai="020.bam/{sample}.bam.bai"
    shell:
        """
        test -s {input.bam}
        samtools index {input.bam}

        echo "testing BAM index..."
        test -s {output.bai}
        """
        
rule samtools_flagstat:
    input:
        bam="020.bam/{sample}.bam"
    output:
        stats="020.bam/{sample}.flagstat"
    shell:
        """
        test -s {input.bam}
        samtools flagstat {input.bam} > {output.stats}

        echo "testing flagstat..."
        test -s {output.stats}
        """

rule bcftools_call:
    input:
        bams=expand("020.bam/{sample}.bam", sample=SAMPLES),
        bais=expand("020.bam/{sample}.bam.bai", sample=SAMPLES)
    output:
        vcf="030.vcf/raw_snps.vcf"
    shell:
        """
        for bam in {input.bams}; do test -s "$bam"; done

        mkdir -p 030.vcf
        bcftools mpileup -Ou -f {DB} {input.bams} | bcftools call -mv -Ov -o {output.vcf}
        
        echo "testing bcftools_call..."
        grep -q "#CHROM" {output.vcf}
        """

rule vt_clean:
    input:
        vcf="030.vcf/raw_snps.vcf"
    output:
        vcf="040.vcf/clean_snps.vcf"
    shell:
        """
        test -s {input.vcf}
        mkdir -p 040.vcf
        cat {input.vcf} \
            | vt decompose - \
            | vt normalize -n -r {DB} - \
            | vt uniq - \
            | vt view -f "QUAL>20" -h - \
            > {output.vcf}

        echo "testing vt_clean..."
        test -s {output.vcf}
        """

rule snpeff_annotate:
    input:
        vcf="040.vcf/clean_snps.vcf"
    output:
        vcf="100.final/snps.annotated.vcf"
    shell:
        """
        test -s {input.vcf}
        mkdir -p 100.final
        java -Xmx3400m -jar {SNPEFF_JAR} eff hg38 \
            -dataDir {SNPEFF_DB} \
            {input.vcf} > {output.vcf}
        
        echo "testing snpeff_annotate..."

        test -s {output.vcf}
        grep -q -m 1 -v "^#" {output.vcf}
        grep -q -m 1 "ANN=" {output.vcf}
        """

rule vcf_to_tsv:
    input:
        vcf="100.final/snps.annotated.vcf"
    output:
        tsv="100.final/snps.annotated.tsv",
        error_log="100.final/snps.errors_summary.txt"
    shell:
        """
        test -s {input.vcf}

        java -jar SnpSift.jar extractFields {input.vcf} CHROM POS REF ALT ANN[*].EFFECT ANN[*].BIOTYPE ANN[*].IMPACT ANN[*].ERRORS > {output.tsv}
        
        echo "testing vcf_to_tsv..."

        test -s {output.tsv}
        head -n 1 {output.tsv} | grep -q "^CHROM"
        awk -F'\\t' 'NR>1 && $8 != ""' {output.tsv} > {output.error_log} || true
        ERR_COUNT=$(wc -l < {output.error_log})
        echo "found $ERR_COUNT variants with annotation errors/warnings"
        """
