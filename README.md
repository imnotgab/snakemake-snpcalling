# **Variant Calling and Annotation Pipeline**
This project contains an automated pipeline based on the Snakemake workflow management system, designed for processing Next-Generation Sequencing (NGS) data. 
It performs an end-to-end analysis from quality control of raw reads, through mapping to the reference genome (hg38 / chr9), to variant calling, cleaning, and annotation of genetic variants (SNPs/Indels).

# **Dependencies**
To run this pipeline, you need the following tools installed (preferably via a Conda environment):

- Snakemake (workflow management)
  
- FastQC (quality control)
  
- BWA (read mapping)
  
- Samtools (BAM file operations)
  
- Bcftools (variant calling)
  
- vt (variant decomposition, normalization, and filtering)
  
- Java (required for Java-based tools)
  
- snpEff (variant annotation)
  
- SnpSift (field extraction from VCF files)

# **Project Structure & Input Data**
The script expects the raw FASTQ read files to be located in the same directory (e.g., TLE66_N.fastq and TLE66_T.fastq).

## Main configuration variables in the script:

- **SAMPLES**: List of analyzed samples (currently: TLE66_N, TLE66_T).
  
- **DB**: Path to the indexed reference genome (currently set to chromosome 9: hg38 / chr9.fa).
  
- **SNPEFF_JAR / SNPEFF_DB**: Paths to the executables and databases for the snpEff tool.
  

# **Pipeline Steps**
The pipeline consists of the following steps:

- **fastq**c: Performs quality control on raw FASTQ files and generates HTML reports.
  
- **bwa_mem_and_sort**: Maps reads to the reference genome using the BWA-MEM algorithm and pipes them to Samtools for sorting and BAM format conversion.

- **samtools_index & samtools_flagstat**: Creates index files (.bai) and generates basic mapping statistics.

- **bcftools_call**: Combines all samples and identifies raw genetic variants (generates a VCF file).

- **vt_clean**: Cleans the data using the _vt toolset_ and decomposes multiallelic variants, normalizes them, removes duplicates, and filters out variants with a quality score below or equal to 20 (QUAL > 20).

- **snpeff_annotate**: Adds functional annotations (e.g., protein impact, transcript biotype) using the hg38 database.

- **vcf_to_tsv**: Extracts selected columns (CHROM, POS, REF, ALT, EFFECT, BIOTYPE, IMPACT, ERRORS) from the annotated VCF file and saves them in a readable tabular format (TSV). Additionally, it filters potential annotation errors into a separate log file without interrupting the analysis.

# **Outputs**
The script will generate the following files:

- **010.fastqc/** - Quality control reports (.html, .zip).

- **020.bam/** - Sorted and indexed BAM files and .flagstat statistics.

- **030.vcf/** - Raw variants in VCF format.

- **040.vcf/** - Cleaned and filtered variants.

- **100.final/** - Final results:

  - **snps.annotated.tsv** - A table ready for analysis in R/Python.

  - **snps.errors_summary.txt** - A log containing warnings and errors generated during the annotation step.
 
# **Usage**
To run the full workflow in standard mode:

```
Bash
#Check what tasks will be executed without actually running them:
snakemake -n

#Actual run (e.g., using 4 CPU cores or more if needed):
snakemake --cores 4
```
