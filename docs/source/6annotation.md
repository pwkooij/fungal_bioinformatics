# Genome annotation
As mentioned during the whole workshop, also here there are many options, however, in the last few years a group of smart people put a nice package together specifically for analyzing fungal genomes named funannotate (https://funannotate.readthedocs.io/).

Genome annotation benefits from RNAseq results but this is not necessary. We can still get good predictions of genes with just a genome assembly.

## Masking repeats
Funannotate has various steps and the first one is to softmask the repeat sections in the genome (it will basically put the repeat sections in the fasta file in small letters and the rest in capital letters):
```
funannotate mask -i PATH/TO/ASSEMBLY --cpus 12 -o PATH/TO/MASKED_ASSEMBLY.fasta
```

## Gene prediction
Now the repeats are masked we can predict genes in the assembly. Funannotate will use AUGUSTUS for this which has pre-trained models for a range of species. We will have to select the species closest to the species we are analyzing. To check the species available we can call up a table with `funannotate species`
Once we have selected the species closest to our sample we can perform the prediction, e.g., `coprinus`:
```
funannotate predict -i PATH/TO/MASKED_ASSEMBLY.fasta -o PATH/TO/OUTPUT --species "YOUR SPECIES NAME" --strain OPTIONAL --busco_seed_species coprinus --cpus 12
```

## Gene identification
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