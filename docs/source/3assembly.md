# 3. Genome assembly

## Illumina only datasets
There are many many assemblers, but for (basidiomycetic) fungi I recommend ABySS:

```
abyss-pe np=8 k=31 B=2G name=path/to/output/prefix in='path/to/reads_1.fastq.gz path/to/reads_2.fastq.gz' &> path/to/verbose.log
```

> **_NOTE:_**
>Since it is always good to test your data, it is recommended to test different kmer sizes and pick the best one.
>The best way to do this is by write a loop script. For example, if you would like to test 15 different kmers with steps of 10 from 20 to 160:
>```
>declare -a k_list=("20" "30" "40" "50" "60" "70" "80" "90" "100" "110" "120" "130" "140" "150" "160")
>
>for k in "${k_list[@]}"; do
>    mkdir PATH/TO/OUTPUT_DIRECTORY/k$k
>    abyss-pe np=8 k=$k -C PATH/TO/OUTPUT_DIRECTORY/k$k name=Ccos_symbiont_k$k in='PATH/TO/forward_paired_trimmed.fastq.gz PATH/TO/reverse_paired_trimmed.fastq.gz' se='PATH/TO/unpaired_trimmed.fastq.gz'
>done
>```

In case you would like to try SPAdes:
```
spades.py --pe1-1 path/to/lib1/reads_1.fastq.gz --pe1-2 path/to/lib1/reads_2.fastq.gz --pe2-1 path/to/lib2/reads_1.fastq.gz --pe2-2 path/to/lib2/reads_2.fastq.gz -o path/to/output -t XX -m XX
```

## Oxford Nanopore only datasets
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

## Hybrid assembly
Sometimes you are lucky and have both Illumina and long reads! Let's use MaSuRCA for this:
```
masurca -t XX -i PATH/TO/illumina_1.fastq.gz,PATH/TO/illumina_2.fastq.gz -r PATH/TO/ont.fastq -p PATH/TO/OUTPUT/ 2>&1 | tee log_masurca.txt
```
