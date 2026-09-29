## Readme for phylogenetic analyses

* `bombina_typeI.fasta` and `bombina_typeII.fasta` contain keratins obtained from HMM search.

`type1_keratins_modified.fasta` and `type2_keratins_modified.fasta` are intermediate filament rod domain alignments produced by adding the predicted proteins for every type I and type II keratin identified in the Bombina bombina genome and the keratins previously identified in Leptobrachium leishanense genome (Li et al. 2019) to the vertebrate keratin type I and type II intermediate filament rod domain alignments of Carron et al. (2024) using the `--keeplength` option in MAFFT (v. 7.526). 

Maximum-likelihood phylogenies were inferred separately for the type I and type II keratin datasets in IQ-TREE (v. 3.1.1). ModelFinder was used to select the best-fitting amino acid substitution model according to the Bayesian information criterion (BIC). Node support was assessed using 1,000 ultrafast bootstrap (UFBoot) replicates with bootstrap NNI optimization.

```
iqtree3 \
  -s type1_keratins_modified.fasta \
  -m MFP \
  -mset LG,WAG,JTT,VT,Dayhoff,Blosum62 \
  -mrate E,G,I,I+G \
  -mfreq FU,F \
  -bb 1000 \
  -bnni \
  -alrt 1000 \
  -nt AUTO

iqtree3 \
  -s type2_keratins_modified.fasta \
  -m MFP \
  -mset LG,WAG,JTT,VT,Dayhoff,Blosum62 \
  -mrate E,G,I,I+G \
  -mfreq FU,F \
  -bb 1000 \
  -bnni \
  -alrt 1000 \
  -nt AUTO
```

IQtree output tree files were visualized using phylogeny_plots.R
