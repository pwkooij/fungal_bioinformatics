# 1. Data Quality Check (QC)

Because sequence reads are generated differently by different sequencers we will also have to use different tools to analyse the data quality for these.

## Illumina QC
For Illumina (and I think PacBIO - CHECK!) we can use, e.g., FastQC:
```
fastqc -t XX -o fastqc file_1.fq.gz file_2.fq.gz
``` 

## Oxford Nanopore QC
For ONT data we can use the R based software MinIONQC
```
MinIONQC.R -i sequencing_summary.txt -o MinIONQC_results/ -p 10 -f 'pdf'
```
