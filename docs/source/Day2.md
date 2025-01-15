# Day 2: Wednesday January 22

## Program for today

- Genome Assembly QC
- Polishing
- Completeness check
- Genome annotation

Hopefully, our assemblies finished during the night, so now it is time to see what we got, how the quality is, polish the assembly, and do some gene annotations!
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
>      - fungi_odb10
         - ascomycota_odb10
             - dothideomycetes_odb10
                 - capnodiales_odb10
                 - pleosporales_odb10
             - eurotiomycetes_odb10
                 - chaetothyriales_odb10
                 - eurotiales_odb10
                 - onygenales_odb10
             - leotiomycetes_odb10
                 - helotiales_odb10
             - saccharomycetes_odb10
             - sordariomycetes_odb10
                 - glomerellales_odb10
                 - hypocreales_odb10
         - basidiomycota_odb10
             - agaricomycetes_odb10
                 - agaricales_odb10
                 - boletales_odb10
                 - polyporales_odb10
             - tremellomycetes_odb10
         - microsporidia_odb10
         - mucoromycota_odb10
             - mucorales_odb10

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
