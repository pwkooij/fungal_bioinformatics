# Day 1: Tuesday January 21

## Program for today

- Brief introduction to polyphasic taxonomy
- Data Quality Check
- Data Trimming
- Genome Assembly

Will use a lot of Linux command line, for reminders have a look at the Help page for technical stuff.

References to the papers about the software and other helpful papers you can find in the References chapter.
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

## IF NECESSARY: basecalling of ONT reads
For Illumina and PacBio, the sequences arrive already basecalled by the machines, i.e., fluorescence information is already translated into DNA. However, for ONT data this is not always the case, especially when sequencing yourself. In this case, we need to translate the electrical current information, or squiggles, that arrive in FAST5 files, into DNA/FASTQ. For newer flowcells and kits we can use either the older software `guppy` or the newer (and better) `dorado`. For older flowcells and kits we have to resort to `guppy`.

To run either we will have to know what flowcell we used, what kit we used, and what settings for the sequencer we used, in order to select the correct model for basecalling, e.g., `dna_r9.5_450bps.cfg` -> DNA sequences with the r9.5 flowcell.

Guppy:
```
guppy_basecaller -i PATH/TO/FAST5_FILES/ -s PATH/TO/OUTPUT/FASTQ/ -c dna_r9.5_450bps.cfg --device cuda:0 --compress_fastq -r
```

Using `dorado` is a bit more elaborate, we first have to put all the FAST5 files together and create a POD5 file from the grouped FAST5 files.

```
single_to_multi_fast5 -i PATH/TO/FAST5_FILES/ -s PATH/TO/CLUSTERED_FAST5/ -t 10 --recursive

pod5 convert fast5 PATH/TO/CLUSTERED_FAST5/ --output PATH/TO/POD5/ -t 10
```

Next, we can start up an analysis server for `dorado`, selecting the correct model, e.g., `dna_r9.4.1_450bps_hac.cfg`-> DNA, on a R9.4 flowcell:
```
dorado_basecall_server --log_path PATH/TO/dorado_server_logs/ --config dna_r9.4.1_450bps_hac.cfg -p auto --device cuda:0
```
This will give us a temporary URL for the port we created which we need for the final step, e.g., `ipc:///tmp/704d-df70-c964-8689`.

And finally, we can do the basecalling of the data:
```
# check port address
ont_basecall_client --config dna_r9.4.1_450bps_hac.cfg --input_path PATH/TO/POD5/ --save_path PATH/TO/DORADO_OUTPUR/ --port ipc:///tmp/704d-df70-c964-8689
```
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
abyss-pe np=8 k=31 B=2G name=path/to/output/prefix in='path/to/reads_1.fastq.gz path/to/reads_2.fastq.gz' &> path/to/verbose.log
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
>```
>smartdenovo -h
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
