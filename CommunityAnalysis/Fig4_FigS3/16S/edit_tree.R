#install.packages("ape")
library("ape")
library("gplots")
library("ggplot2")
library("ggtree")
#library("phangorn")
#library("gridExtra")
#library("dada2")
#library("msa")
#library("phyloseq")
#source("https://bioconductor.org/biocLite.R")
#biocLite('msa')
#library("seqinr")
#biocLite('DECIPHER')
#library("DECIPHER")

#Before this can be run, the sequences from the file made by the 'get_abundant_sequences.py' script should be aligned to the SILVA database using the ARB aligner:
#https://www.arb-silva.de/aligner/
#The results of this alignment then need to be turned into a tree using QIIME with the command:
#make_phylogeny.py -i aligned_fasta_from_ARB.fasta -o ARB_aligned_mid_rooted.dnd' -r midpoint
#By default QIIME uses fasttree to construct the tree

tree <- read.tree("/Users/u1560915/Documents/OneDrive/PhD_Plastic_Oceans/Experiments/MiSeq_Dada/tree/16S_ARB_aligner/ARB_aligned_mid_rooted.dnd")
pdf(filename="/Users/u1560915/Documents/OneDrive/PhD_Plastic_Oceans/Experiments/MiSeq_Dada/tree/16S_tree_test/new_tree_align.pdf")
plt <- ggtree(tree)
plt + geom_tiplab(size=3, color="purple")+ geom_text2(aes(label=node))
plt <- flip(plt, 71, 13)
plt + geom_tiplab(size=3, align=TRUE, linesize=.5)
#+geom_treescale(0,-1)

tree$tip.label

dev.off()