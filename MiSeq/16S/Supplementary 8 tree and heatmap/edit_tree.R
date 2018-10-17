install.packages("ape")
library("ape")
library("gplots")
library("ggplot2")
library("ggtree")
library("treeio")
vignette("Importer", package="treeio")
tree <- read.tree("/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Supplementary 8 tree and heatmap/16S_new_fasta_0.5%.dnd")
pdf("/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/16S/Supplementary 8 tree and heatmap/new_tree_no_length_colors.pdf")

plt <- ggtree(tree, branch.length="none")
#, branch.length="none"
#plt + geom_tiplab(size=3, color="purple")
plt
plt + geom_hilight(node=74, fill="darkorange1") +  geom_hilight(node=70, fill="darkorange1") + geom_hilight(node=66, fill="blueviolet") + geom_hilight(node=65, fill="cadetblue1") + geom_hilight(node=68, fill="deeppink2") + geom_hilight(node=60, fill="gold1") + geom_hilight(node=46, fill="blue1") + geom_hilight(node=77, fill="green3") + geom_hilight(node=86, fill="gold1")
#geom_hilight(node=72, fill="steelblue", alpha=.6) +
#geom_hilight(node=17, fill="steelblue", alpha=.6) +
#+ geom_text2(aes(subset=!isTip, label=node)) 
plt

dev.off()
print(col2rgb("green3"))