# 4. Genome Assembly QC
We can check our assemblies in a variety of ways, e.g., compare it to reference genomes, check the overall statistics, see if there is any contamination that needs to be removed, etc

## Assembly statistics
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
## Completeness check (before polishing)
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

## Contamination and overall statistics check
Even if we try to work as clean as possible we sometimes still pick up contaminations that end up in our sequence results. Luckily there are tools that can help identify contaminated reads and select to remove them. Blobtools is one such tool: it uses BLAST, Diamond, minimap2, and BUSCO to give overall statistics and matches to the closest organisms for each of the reads.
The software is quite extensive, so we might skip this but I will leave the codes here for you to use (if there is time or in the future).

Also, I have never run this on a server and it seems they have some instructions for this: https://blobtoolkit.genomehubs.org/pipeline/pipeline-tutorials/running-the-pipeline-on-a-cluster/

The full package of blobtools, know as blobtoolkit, exists of several steps that can be performed in whichever order you prefer with the exception of setting up the results folder and creating the base file for the analyses:
### Create meta data yaml file
First we will need to create a file called `meta.yaml`which contains some basic information about your sample, e.g.;
```
record_type: scaffold 
taxon: 
	name: Leucoagaricus gongylophorus 
	taxid: 79220 
	phylum: Basidiomycota
```
As previously, you can use your favourite text editor for this (e.g., nano, gedit, vim, etc)
### Create a BlobDir
Next we need to create a results folder called the BlobDir and add the information from the `meta.yaml`file as well as some information regarding the taxonomy of the sample. For this you will need to know the NCBI taxid (https://www.ncbi.nlm.nih.gov/taxonomy):
```
blobtools create --fasta PATH/TO/ASSEMBLY_FILE --meta meta.yaml --taxid XXXX --taxdump PATH/TO/taxdump PATH/TO/OUTPUT/DIRECTORY
```
And we will merge all information together into a single file:
```
cat PATH/TO/BlobDir/meta.json
```
### Use BLAST and UNIPROT to identify contaminations
This part works best when you have a local copy of the BLAST and UNIPROT databases, however, this takes up almost a 1TB...
#### BLAST 
BLAST is a fantastic tool but it isn't very efficient in using its resources, i.e., blasting my assemblies generally takes about 3 days with 10 threads:
```
blastn -db nt -query PATH/TO/ASSEMBLY_FILE -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 10 -out PATH/TO/BlobDir/blast.out
```
#### UNIPROT
UNIPROT is faster than BLAST but looks at proteins only, and, therefore, tends to have less hits:
```
diamond blastx --query PATH/TO/ASSEMBLY_FILE --db PATH/TO/uniprot/reference_proteomes.dmnd --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 10 > PATH/TO/BlobDir/diamond.out
```
#### Adding BLAST and UNIPROT results to BlobDir
Now we have generated the BLAST and UNIPROT results we need to add them to the overall BlobDir results:
```
blobtools add --hits PATH/TO/BlobDir/blast.out --hits PATH/TO/BlobDir/diamond.out --taxrule bestsumorder --taxdump PATH/TO/taxdump PATH/TO/BlobDir
```
### Adding mapping coverage
Next we will mapping coverage information for our assembly by mapping our reads back on to the assembly. This will show use how much coverage we will have for each contig.
#### Using minimap2 to map reads
```
# SHORT READS
minimap2 -ax sr -t 10 PATH/TO/ASSEMBLY_FILE PATH/TO/READS/illumina_1.fastq.gz PATH/TO/READS/illumina_2.fastq.gz | samtools sort -@10 -O BAM -o PATH/TO/BlobDir/NAME.shortreads.bam -

# LONG READS -> in case of PacBio use -ax map-pb
minimap2 -ax map-ont -t 10 PATH/TO/ASSEMBLY_FILE PATH/TO/READS/ont.fastq | samtools sort -@10 -O BAM -o PATH/TO/BlobDir/NAME.longreads.bam -
```
#### Adding mapping coverage of short and/or long reads to BlobDir
And again we will add the results to the BlobDir:
```
# SHORT READS
blobtools add --cov PATH/TO/BlobDir/NAME.shortreads.bam PATH/TO/BlobDir

# LONG READS
blobtools add --cov PATH/TO/BlobDir/NAME.longreads.bam PATH/TO/BlobDir
```
### Add BUSCO scores
#### Run BUSCO if necessary
If BUSCO wasn't done yet, now is the time to do that, otherwise you can skip this step and immediately add your BUSCO results to the BlobDir:
```
busco -i PATH/TO/ASSEMBLY_FILE -o PATH/TO/OUTPUT -l YOUR_LINEAGE -m geno -c XX
```
#### Add to BUSCO results to BlobDir
```
blobtools add --busco PATH/TO/OUTPUT/run_YOUR_LINEAGE/full_table.tsv PATH/TO/BlobDir
```
### Add links to BlobDir
The next step will just add more background information to the results output and is optional:
```
blobtools add --link taxon.taxid.ENA="https://www.ebi.ac.uk/ena/ EBI ENA LINK OF YOUR SPECIES" --link taxon.name.Wikipedia="WIKIPEDIA LINK OF YOUR SPECIES" PATH/TO/BlobDir
```
### Open dataset in BlobToolKit Viewer
Now we have collected all information and it is time to visualize them! This will create a local URL you can open in a browser to look at Blobtools beautiful graphs:
```
blobtools view --remote PATH/TO/BlobDir
```
