# 5. Polishing the assembly
One way of improving the genome assembly is by comparing the obtained contigs to the reads used to create them, which we call polishing. 
Once again, depending on the sequencing technique, we will have to use different software to do this.

## Polishing Illumina and/or hybrid data
For this we can use nextPolish (https://nextpolish.readthedocs.io/en/latest/index.html). To use nextPolish we will have to create a list of the read files:

```
ls PATH/TO/reads1_R1.fq PATH/TO/reads1_R2.fq > PATH/TO/WORKDIR/sgs.fofn
ls PATH/TO/ont_reads.fq > PATH/TO/WORKDIR/lgs.fofn
```

We will also need to create a config file which will tell the program where to find the files, where to put the output, and what options to use:
```
nano PATH/TO/WORKDIR/run.cfg
```

This will open a new file with the name run.cfg. There we will fill out the following information:
```
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = PATH/TO/ASSEMBLY_FILE
genome_size = auto
workdir = PATH/TO/WORKDIR/
polish_options = -p {multithread_jobs}

[sgs_option] # Short-read settings

sgs_fofn = PATH/TO/sgs.fofn
sgs_options = -max_depth 100 -bwa

[lgs_option] # Long-read settings
lgs_fofn = PATH/TO/lgs.fofn
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont # in case of PacBio use map-pb
```
Now we have set up the config file and the lists with reads and it is time to start the analysis!
```
nextPolish run.cfg
```

## Polishing Oxford Nanopore data
Because the error rate in ONT reads is relatively high (compared to other sequencing methods) polishing the assembly needs a bit more work. What is generally recommended is perform Racon 4x and then Medaka. In some cases it is possible to run the software Homopolish after that as well, but this software relies heavily on already available genomes from closely related species (i.e., in my case, working with *Leucocoprineae* this didn't work).

### Polishing four rounds with Racon
To run Racon, we will prepare our data first with minimap2 followed by polishing the assembly with Racon:
```
minimap2 -a -t XX [ASSEMBLY FILE] [READS FILE] > OUTPUT.sam
racon -m 8 -x -6 -g -8 -w 500 -t XX [READS FILE] OUTPUT.sam [ASSEMBLY FILE] > OUTPUT.fasta
```

> **For Racon settings we have to look at the help output:**
>```
>racon {options} [sequences]  [overlaps]  [target sequences] # default output is stdout
>   [sequences] # input file in FASTA/FASTQ format (can be compressed with gzip) containing sequences used for correction
>    [overlaps] # input file in MHAP/PAF/SAM format (can be compressed with gzip) containing overlaps between sequences and target sequences
>    [target sequences] # input file in FASTA/FASTQ format (can be compressed with gzip) containing sequences which will be corrected
>
>options:
>    -u, --include-unpolished # output unpolished target sequences
>    -f, --fragment-correction # perform fragment correction instead of contig polishing (overlaps file should contain dual/self overlaps!)
>    -w, --window-length [int] # default: 500, size of window on which POA is performed
>    -q, --quality-threshold [float] # default: 10.0, threshold for average base quality of windows used in POA
>    -e, --error-threshold [float] # default: 0.3, maximum allowed error rate used for filtering overlaps
>    --no-trimming # disables consensus trimming at window ends
>    -m, --match [int] # default: 3, score for matching bases
>    -x, --mismatch [int] # default: -5, score for mismatching bases
>    -g, --gap [int] # default: -4, gap penalty (must be negative)
>    -t, --threads [int] # default: 1, number of threads
>    --version # prints the version number
>    -h, --help # prints the usage
>
>only available when built with CUDA:
>    -c, --cudapoa-batches [int] # default: 0, number of batches for CUDA accelerated polishing per GPU
>    -b, --cuda-banded-alignment # use banding approximation for polishing on GPU. Only applicable when -c is used.
>    --cudaaligner-batches [int] # default: 0, number of batches for CUDA accelerated alignment per GPU
>    --cudaaligner-band-width [int] # default: 0, Band width for cuda alignment. Must be >= 0. Non-zero allows user defined band width, whereas 0 implies auto band width determination.
>```
   
>**NOTE:**
>Sometimes the assembly file has a different file extension that `minimap2`is not able to recognize, e.g., in the case of SMARTdenovo, where the assembly file ends in `.cns`. In that case we will have to add `.fasta` as the file extension.
>```
>cp PATH/TO/ASSEMBLY.zmo.cns PATH/TO/WORKDIR/
>mv PATH/TO/WORKDIR/ASSEMBLY.zmo.cns PATH/TO/WORKDIR/ASSEMBLY.zmo.cns.fasta
>```

Because it is recommended to perform four rounds of Racon polishing we will have to redo the analysis 3x each with the newly produced assembly produced in the previous round:
```
minimap2 -a -t XX [RACON1 ASSEMBLY FILE] [READS FILE] > OUTPUT.sam
racon -m 8 -x -6 -g -8 -w 500 -t XX [READS FILE] OUTPUT.sam [RACON1 ASSEMBLY FILE] > OUTPUT.fasta
```

### Polishing with Medaka
To run Medaka, it is important to know what Medaka model to use. To figure this out we have to look at what we used to obtain the data.

The model are named with the following scheme:
```
{pore}_{device}_{caller variant}_{caller version}
```
For example, if the pore/flowcell is r941, the device is MinION (min), you used the high-accuracy model (high), and the guppy version was 4.15. Choose the model, that is closest to that basecaller version.

With de medaka_consenus help we find:
```
-m  medaka model, (default: r1041_e82_400bps_sup_v5.0.0).
        Choices: r103_fast_g507 r103_hac_g507 r103_sup_g507 r1041_e82_260bps_fast_g632 r1041_e82_260bps_hac_g632 r1041_e82_260bps_hac_v4.0.0 r1041_e82_260bps_hac_v4.1.0 r1041_e82_260bps_joint_apk_ulk_v5.0.0 r1041_e82_260bps_sup_g632 r1041_e82_260bps_sup_v4.0.0 r1041_e82_260bps_sup_v4.1.0 r1041_e82_400bps_bacterial_methylation r1041_e82_400bps_fast_g615 r1041_e82_400bps_fast_g632 r1041_e82_400bps_hac_g615 r1041_e82_400bps_hac_g632 r1041_e82_400bps_hac_v4.0.0 r1041_e82_400bps_hac_v4.1.0 r1041_e82_400bps_hac_v4.2.0 r1041_e82_400bps_hac_v4.3.0 r1041_e82_400bps_hac_v5.0.0 r1041_e82_400bps_sup_g615 r1041_e82_400bps_sup_v4.0.0 r1041_e82_400bps_sup_v4.1.0 r1041_e82_400bps_sup_v4.2.0 r1041_e82_400bps_sup_v4.3.0 r1041_e82_400bps_sup_v5.0.0 r104_e81_fast_g5015 r104_e81_hac_g5015 r104_e81_sup_g5015 r104_e81_sup_g610 r941_e81_fast_g514 r941_e81_hac_g514 r941_e81_sup_g514 r941_min_fast_g507 r941_min_hac_g507 r941_min_sup_g507 r941_prom_fast_g507 r941_prom_hac_g507 r941_prom_sup_g507 r941_sup_plant_g610
```

So in this case we need to choose `r941_min_hac_g507`.

Now we are ready to run Medaka:
```
medaka_consensus -i PATH/TO/ONT.fastq.gz -d PATH/TO/RACON4_ASSEMBLY.fasta -o PATH/TO/OUTPUT -t XX -m MODEL
```

### Polishing with Homopolish

INFO TO COME