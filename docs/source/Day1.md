# Day 1: Tuesday January 21

## Program for today

- Brief introduction to polyphasic taxonomy
- Data Quality Check
- Data Trimming
- Genome Assembly

## Data Quality Check (QC)

Because sequence reads are generated differently by different sequencers we will also have to use different tools to analyse the data quality for these.

### Illumina QC
For Illumina (and I think PacBIO - CHECK!) we can use, e.g., FastQC:
```
fastqc -t XX -o fastqc file_1.fq.gz file_2.fq.gz
``` 

### Oxford Nanopore QC
For ONT data we can use the R based software MinIONQC
```
MinIONQC.R -i sequencing_summary.txt -o MinIONQC_results/ -p 10 -f 'pdf'
```

## Data Trimming

### Illumina data trimming
Also here, there are many options, e.g., Trimmomatic, fastp, and more
fastp:
```
fastp --in1 file_1.fq.gz --in2 file_2.fq.gz --out1 illumina_1.fastq.gz --out2 illumina_2.fastq.gz --unpaired1 illumina_u.fastq.gz --unpaired2 illumina_u.fastq.gz
``` 
Trimmomatic:
```
java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz
```

### Oxford Nanopore data trimming
For this we will use here porechop. Sometimes it is necessary to join several fastq files into one:
```
cat fastq_runid_*.fastq.gz > porechop/all.fastq.gz
cd porechop
porechop -i all.fastq.gz -o porechop_all.fastq.gz --threads 10 > porechop.out
```