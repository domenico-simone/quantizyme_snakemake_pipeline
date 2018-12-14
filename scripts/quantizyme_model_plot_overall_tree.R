required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs)
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))
suppressMessages(library(phangorn))
suppressMessages(library(ape))
suppressMessages(library(optparse))

option_list = list(
    make_option("--tree", type="character", default="hc_tree_2.pdf",
            help="Tree file version 2 - pdf"),
    make_option("--tree_plot", type="character", default="tree2.new",
            help="Tree file version 2 - txt"),
    make_option(c("-p", "--projectID"), type="character", default=NULL,
            help="input fasta file [default= %default] (MANDATORY)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

overall_tree = opt$tree
tree_plot = opt$tree_plot
project.id = opt$projectID

cluster2 = read.tree(overall_tree)
if(class(cluster2) != "phylo") stop(" Could not calculate phylogenetic tree by neighbor-joining.\n\n")

if(cluster2$Nnode < 25) {
    labels = TRUE
} else {
    labels = FALSE
}

#X11(width=6, height=6)
jpeg(tree_plot)
plot.phylo(cluster2, type = "unrooted", show.tip.label = labels, main = paste("Phylogenetic tree, project", project.id), font.main = 1)
invisible(dev.off())
# fn = paste(project.id, "hc_tree_2.png", sep="_")
# outpdf7 = file.path(getwd(), outfolder, fn)
# invisible(dev.copy2pdf(file = outpdf7, out.type = "pdf")); invisible(dev.off())
cat(paste("  Phylogenetic tree (version 2) saved to", tree_plot, "\n"))
##### link plot
# insert.text(outhtml, text="Phylogenetic tree, radial format (nj, ape):", color="black")
# set.link(outhtml, outpdf7, caption="pdf")
# insert.blanks(outhtml, 2)
# link tree in Newick format:
# fn = "tree2.new"
# tree.file = file.path(outdir, outtxt7)
# #tree.file = file.path(getwd(), outfolder, fn)
# write.tree(phy = cluster2, file = tree.file)
# cat(paste("  Neighbor joining tree (version 2) saved to '", outtxt7, "'\n\n"))
