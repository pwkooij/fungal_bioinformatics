# Day 2: Wednesday January 22

## Program for today

- Genome Assembly QC
- Polishing
- Completeness check
- Genome annotation

Hopefully, our assemblies from yesterday finished during the night, so now it is time to see what we got, how the quality is, polish the assembly, and do some gene annotations!
## Genome Assembly QC
We can check our assemblies in a variety of ways, e.g., compare it to reference genomes, check the overall statistics, see if there is any contamination that needs to be removed, etc

### Assembly statistics
We will start with some simple statistics such as N50, number of contigs, contig size, etc.
For this we use a simple perl script developed by Dr. Heath O'Brien (https://github.com/hobrien)
```
ContigStats -i ASSEMBLY_FILE -f FILE_FORMAT > output.log
```

When we then open the `output.log` file with a text editor (e.g., nano, gedit, vim), we will see something similar to this:
```
N50: 126,880  
Number contigs: 693  
Max_contig: 571,286  
Min_contig: 869  
Total_length: 37,632,407  
Average_contig: 54,304  
Num_gaps: 80  
Total_gap_length: 323
```

These are some basic statistics, but if we want to get more details or compare it to a reference genome of the same or closely related species, we will have to use the tool called QUAST - QUality ASsessment Tool:

```
quast [SEQUENCE FILE] -r [REFERENCE ASSEMBLY optional] -o [OUTPUT_NAME] -t XX --nanopore [ONT READS FILE] -1 [ILLUMINA 1 fastq.gz] -2 [ILLUMINA 2 fastq.gz]
```

This gives a HTML output that we can open in a browser with various statistics.
### Completeness check (before polishing)
Statistics on size and number of fragments can give a good idea of the quality, but don't tell everything. You might have beautiful statistics, but if the DNA of the fragments is of low quality, it is still a bad assembly! So how do we check this? But comparing the contents to fixed databases.

For this we use BUSCO, a tool that looks in a database of single copy genes for a particular group of organisms, and checks if it can find those genes in your assembly and if they are presented one time or are duplicated or fragmented

```
busco -i [ASSEMBLY FILE] -l [DATASET NAME] -o [OUTPUT_NAME] -m geno -c 10
```

> **List of fungal datasets available at the moment**:
>- fungi_odb10
>- ascomycota_odb10
>- dothideomycetes_odb10
>- capnodiales_odb10
>- pleosporales_odb10
>- eurotiomycetes_odb10
>- chaetothyriales_odb10
>- eurotiales_odb10
>- onygenales_odb10
>- leotiomycetes_odb10
>- helotiales_odb10
>- saccharomycetes_odb10
>- sordariomycetes_odb10
>- glomerellales_odb10
>- hypocreales_odb10
>- basidiomycota_odb10
>- agaricomycetes_odb10
>- agaricales_odb10
>- boletales_odb10
>- polyporales_odb10
>- tremellomycetes_odb10
>- microsporidia_odb10
>- mucoromycota_odb10
>- mucorales_odb10

Example of an output:
```
C:95.3% (S:94.7%, D:0.6%), F:0.7%, M:4.0%, n:3870, E:5.7%       
3689    Complete BUSCOs (C)     (of which 209 contain internal stop codons)
3666    Complete and single-copy BUSCOs (S)       
23      Complete and duplicated BUSCOs (D)        
28      Fragmented BUSCOs (F)                     
153     Missing BUSCOs (M)                        
3870    Total BUSCO groups searched
```

### Contamination and overall statistics check
Even if we try to work as clean as possible we sometimes still pick up contaminations that end up in our sequence results. Luckily there are tools that can help identify contaminated reads and select to remove them. Blobtools is one such tool: it uses BLAST, Diamond, minimap2, and BUSCO to give overall statistics and matches to the closest organisms for each of the reads.
The software is quite extensive, so we might skip this but I will leave the codes here for you to use (if there is time or in the future).

Also, I have never run this on a server and it seems they have some instructions for this: https://blobtoolkit.genomehubs.org/pipeline/pipeline-tutorials/running-the-pipeline-on-a-cluster/

The full package of blobtools, know as blobtoolkit, exists of several steps that can be performed in whichever order you prefer with the exception of setting up the results folder and creating the base file for the analyses:
#### Create meta data yaml file
First we will need to create a file called `meta.yaml`which contains some basic information about your sample, e.g.;
```
record_type: scaffold 
taxon: 
	name: Leucoagaricus gongylophorus 
	taxid: 79220 
	phylum: Basidiomycota
```
As previously, you can use your favourite text editor for this (e.g., nano, gedit, vim, etc)
#### Create a BlobDir
Next we need to create a results folder called the BlobDir and add the information from the `meta.yaml`file as well as some information regarding the taxonomy of the sample. For this you will need to know the NCBI taxid (https://www.ncbi.nlm.nih.gov/taxonomy):
```
blobtools create --fasta PATH/TO/ASSEMBLY_FILE --meta meta.yaml --taxid XXXX --taxdump PATH/TO/taxdump PATH/TO/OUTPUT/DIRECTORY
```
And we will merge all information together into a single file:
```
cat PATH/TO/BlobDir/meta.json
```
#### Use BLAST and UNIPROT to identify contaminations
This part works best when you have a local copy of the BLAST and UNIPROT databases, however, this takes up almost a 1TB...
##### BLAST 
BLAST is a fantastic tool but it isn't very efficient in using its resources, i.e., blasting my assemblies generally takes about 3 days with 10 threads:
```
blastn -db nt -query PATH/TO/ASSEMBLY_FILE -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 10 -out PATH/TO/BlobDir/blast.out
```
##### UNIPROT
UNIPROT is faster than BLAST but looks at proteins only, and, therefore, tends to have less hits:
```
diamond blastx --query PATH/TO/ASSEMBLY_FILE --db PATH/TO/uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 10 > PATH/TO/BlobDir/diamond.out
```
##### Adding BLAST and UNIPROT results to BlobDir
Now we have generated the BLAST and UNIPROT results we need to add them to the overall BlobDir results:
```
blobtools add --hits PATH/TO/BlobDir/blast.out --hits PATH/TO/BlobDir/diamond.out --taxrule bestsumorder --taxdump PATH/TO/taxdump PATH/TO/BlobDir
```
#### Adding mapping coverage
Next we will mapping coverage information for our assembly by mapping our reads back on to the assembly. This will show use how much coverage we will have for each contig.
##### Using minimap2 to map reads
```
# SHORT READS
minimap2 -ax sr -t 10 PATH/TO/ASSEMBLY_FILE PATH/TO/READS/illumina_1.fastq.gz PATH/TO/READS/illumina_2.fastq.gz | samtools sort -@10 -O BAM -o PATH/TO/BlobDir/NAME.shortreads.bam -

# LONG READS -> in case of PacBio use -ax map-pb
minimap2 -ax map-ont -t 10 PATH/TO/ASSEMBLY_FILE PATH/TO/READS/ont.fastq | samtools sort -@10 -O BAM -o PATH/TO/BlobDir/NAME.longreads.bam -
```
##### Adding mapping coverage of short and/or long reads to BlobDir
And again we will add the results to the BlobDir:
```
# SHORT READS
blobtools add --cov PATH/TO/BlobDir/NAME.shortreads.bam PATH/TO/BlobDir

# LONG READS
blobtools add --cov PATH/TO/BlobDir/NAME.longreads.bam PATH/TO/BlobDir
```
#### Add BUSCO scores
##### Run BUSCO if necessary
If BUSCO wasn't done yet, now is the time to do that, otherwise you can skip this step and immediately add your BUSCO results to the BlobDir:
```
busco -i PATH/TO/ASSEMBLY_FILE -o PATH/TO/OUTPUT -l YOUR_LINEAGE -m geno -c XX
```
##### Add to BUSCO results to BlobDir
```
blobtools add --busco PATH/TO/OUTPUT/run_YOUR_LINEAGE/full_table.tsv PATH/TO/BlobDir
```
#### Add links to BlobDir
The next step will just add more background information to the results output and is optional:
```
blobtools add --link taxon.taxid.ENA="https://www.ebi.ac.uk/ena/ EBI ENA LINK OF YOUR SPECIES" --link taxon.name.Wikipedia="WIKIPEDIA LINK OF YOUR SPECIES" PATH/TO/BlobDir
```
#### Open dataset in BlobToolKit Viewer
Now we have collected all information and it is time to visualize them! This will create a local URL you can open in a browser to look at Blobtools beautiful graphs:
```
blobtools view --remote PATH/TO/BlobDir
```

## Polishing the assembly
One way of improving the genome assembly is by comparing the obtained contigs to the reads used to create them, which we call polishing. 
Once again, depending on the sequencing technique, we will have to use different software to do this.

### Polishing Illumina and/or hybrid data
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

### Polishing Oxford Nanopore data
Because the error rate in ONT reads is relatively high (compared to other sequencing methods) polishing the assembly needs a bit more work. What is generally recommended is perform Racon 4x and then Medaka. In some cases it is possible to run the software Homopolish after that as well, but this software relies heavily on already available genomes from closely related species (i.e., in my case, working with *Leucocoprineae* this didn't work).

#### Polishing four rounds with Racon
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

#### Polishing with Medaka
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

## Genome annotation
As mentioned during the whole workshop, also here there are many options, however, in the last few years a group of smart people put a nice package together specifically for analyzing fungal genomes named funannotate (https://funannotate.readthedocs.io/).

Genome annotation benefits from RNAseq results but this is not necessary. We can still get good predictions of genes with just a genome assembly.

Funannotate has various steps and the first one is to softmask the repeat sections in the genome (it will basically put the repeat sections in the fasta file in small letters and the rest in capital letters):
```
funannotate mask -i PATH/TO/ASSEMBLY --cpus 12 -o PATH/TO/MASKED_ASSEMBLY.fasta
```

Now the repeats are masked we can predict genes in the assembly. Funannotate will use AUGUSTUS for this which has pre-trained models for a range of species. We will have to select the species closest to the species we are analyzing. To check the species available we can call up a table with `funannotate species`
Once we have selected the species closest to our sample we can perform the prediction, e.g., `coprinus`:
```
funannotate predict -i PATH/TO/MASKED_ASSEMBLY.fasta -o PATH/TO/OUTPUT --species "YOUR SPECIES NAME" --strain OPTIONAL --busco_seed_species coprinus --cpus 12
```

This will create a set of gene models from our assembly that we can use to compare to various databases. Using `funannotate` we will compare the gene models with `interproscan`, `eggnog`, `antismash`, and `phobius`. We will start with `interproscan`:
```
#run using docker
funannotate iprscan -i PATH/TO/PREDICT_OUTPUT/predict_results/ -m docker --cpus 12

#run locally (Linux only)
funannotate iprscan -i PATH/TO/PREDICT_OUTPUT/predict_results/ -m local --iprscan_path /my/path/to/interproscan.sh -o OUTPUT
```

Next we will analyse the data with `antismash` and `phobius`. If they are not installed locally, we can do both at the same time remotely:
```
funannotate remote -i PATH/TO/PREDICT_OUTPUT -m phobius antismash -e your-email@domain.edu
```

>**NOTE:**
>This doesn't work properly and is quite slow. We should have both `phobius` and `antismash` installed on the server, so we can run it locally.


>**NOTE:**
>THE BELOW PART IT STILL IN PROGRESS!

For `phobius` locally:
```
phobius -short -i PATH/TO/annotate_misc/genome.proteins.fasta -o annotate_misc/phobius.results.txt -l logfiles/phobius.log 
```

For `antismash` locally:
```
antismash -t fungi -c 64 --databases PATH/TO/ANTISMASH/DATABASES --output-dir PATH/TO/OUTPUT --output-basename PREFIX -i PATH/TO/PREDICT_OUTPUT.gbk
```

Finally, we can combine all the analyses and do the actual annotation. If `eggnog` is installed locally, it will automatically add this analysis as well:
```
funannotate annotate -i PATH/TO/PREDICT_OUTPUT --cpus 64
```

Alternatively, when running the `antismash` and `phobius` locally:
```
funannotate annotate -i PATH/TO/PREDICT_OUTPUT --antismash PATH/TO/ANTISMASH_RESULTS --cpus 64
```