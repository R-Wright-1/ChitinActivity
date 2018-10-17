library("ape")
library("gplots")
library("ggplot2")
library("ggtree")
library("treeio")
vignette("Importer", package="treeio")

tree <- read.tree("/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Supplementary 9 tree and heatmap/18S_new_fasta_0.5%.dnd")
pdf("/Users/u1560915/Documents/GitHub/ChitinActivity/MiSeq/18S/Supplementary 9 tree and heatmap/new_tree_label_color.pdf")
plt <- ggtree(tree, branch.length="none")

plt + geom_hilight(node=171, fill="darkorange1") + geom_hilight(node=175, fill="cadetblue1") + geom_hilight(node=177, fill="blueviolet") + geom_hilight(node=180, fill="cadetblue1") + geom_hilight(node=163, fill="gold1") + geom_hilight(node=159, fill="blue1") + geom_hilight(node=157, fill="green3") + geom_hilight(node=150, fill="red1") + geom_hilight(node=117, fill="darkblue") + geom_hilight(node=113, fill="cornsilk4") + geom_hilight(node=139, fill="darkolivegreen1") + geom_hilight(node=96, fill="darkmagenta") + geom_hilight(node=103, fill="darkred")

#geom_hilight(node=72, fill="steelblue", alpha=.6) +
#geom_hilight(node=17, fill="steelblue", alpha=.6) +
#+ geom_text2(aes(subset=!isTip, label=node)) +  

#, branch.length="none"
#plt + geom_tiplab(size=3, color="purple")
plt
dev.off()

print(col2rgb("darkred"))