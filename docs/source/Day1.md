# Day 1: Tuesday January 21

## Program for today

- Brief introduction to polyphasic taxonomy
- Data Quality Check
- Data Trimming
- Genome Assembly

Will use a lot of Linux command line, for reminders have a look out the [[Tech help|tech section]].

References to the papers about the software and other helpful papers you can find in the [[References]] chapter.
## Brief introduction to polyphasic taxonomy
Some helpful links:
> **To check name status**
> - Species Fungorum: https://www.speciesfungorum.org/
> - Index Fungorum: https://www.indexfungorum.org/
> - Mycobank: https://www.mycobank.org/
> - Fungal Names: https://nmdc.cn/fungalnames/allfungal
> 
> **Other helpful sites**
> - Interactive identification key for seed fungal pathogens: https://keys.lucidcentral.org/keys/v3/seed_fungi/veg_seed_fungi.html
> - Unite: https://unite.ut.ee/
> - Mycoportal: https://www.mycoportal.org/portal/
> - Mould identification: a virtual self assessment: https://www.adelaide.edu.au/mycology/mould-identification-a-virtual-self-assessment

.. note::

   Testing the note this way

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

## Genome assembly

### Illumina only datasets
There are many many assemblers, but for (basidiomycetic) fungi I recommend ABySS:

```
abyss-pe np=8 k=31 name=path/to/output/prefix in='path/to/reads_1.fastq.gz path/to/reads_2.fastq.gz' &> path/to/verbose.log
```

> **_NOTE:_**
>Since it is always good to test your data, it is recommended to test different kmer sizes and pick the best one.

In case you would like to try SPAdes:
```
spades.py --pe1-1 path/to/lib1/reads_1.fastq.gz --pe1-2 path/to/lib1/reads_2.fastq.gz --pe2-1 path/to/lib2/reads_1.fastq.gz --pe2-2 path/to/lib2/reads_2.fastq.gz -o path/to/output -t XX -m XX
```

### Oxford Nanopore only datasets
I have tried various assemblers and so far my favourite is SMARTdenovo:
```
smartdenovo -p PREFIX -c 1 -e zmo -t XX -J 500 PATH/TO/TRIMMED.fastq.gz > FANCYNAME.mak
make -f FANCYNAME.mak
```
> **_NOTE:_**
> SMARTdenovo -h gives:
> ```
> smartdenovo -h
>
>Unknown option: h
>
>Usage: smartdenovo.pl [options] <reads.fa>
>
>Options:
 > -p STR     output prefix [wtasm]
 > -e STR     engine of overlaper, compressed kmer overlapper(zmo), dot matrix overlapper(dmo), [dmo]
>  -t INT     number of threads [8]
 > -k INT     k-mer length for overlapping [16]
  >           for large genome as human, please use 17
 > -J INT     min read length [5000]
 > -c INT     generate consensus, [0]
>```
> Consider if you want to change -J read length to a different size, change the kmer length, or if you want to use their newer dot matrix overlapper

### Hybrid assembly
Sometimes you are lucky and have both Illumina and long reads! Let's use MaSuRCA for this:
```
masurca -t XX -i PATH/TO/illumina_1.fastq.gz,PATH/TO/illumina_2.fastq.gz -r PATH/TO/ont.fastq -p PATH/TO/OUTPUT/ 2>&1 | tee log_masurca.txt
```
