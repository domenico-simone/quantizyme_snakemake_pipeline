required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs)
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))
suppressMessages(library(phangorn))
suppressMessages(library(ape))
suppressMessages(library(optparse))

option_list = list(
    make_option(c("-a", "--alignment_file"), type="character", default=NULL,
        help="alignment file [default= %default] (MANDATORY)", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="input fasta file [default= %default] (MANDATORY)"),
    make_option(c("-n", "--n_subtrees"), type="integer", default=NULL,
        help="number of subtrees [default= %default] (MANDATORY)", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


nrHMM = as.integer(opt$n)
transcript.fasta = opt$i

transcript.fas = read.fasta(transcript.fasta, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)

## +++ SUBTREEing
# insert.line(outhtml, percent=100, size=1)

subtree_by_list = FALSE

if(subtree_by_list) {
	cat(paste("  Subtreeing is conducted according to the list '", subtrees.txt, "'\n"))
	success <- file.copy(subtrees.txt, file.path(getwd(), outfolder), overwrite = TRUE)
	if(!success) stop(paste(" File '", subtrees.txt, "' could not be copied to folder", outfolder, ".\n\n"))
	insert.text(outhtml, text="Subtreeing is conducted according to the list ", color="black")
	set.link(outhtml, file.path(getwd(), outfolder, subtrees.txt), caption=subtrees.txt)
	line.feed(outhtml, 1)

	nrHMM = length(unique(subtrees[,2]))
	subtree = list()
	subtree.fn = list()
	sum = 0
	for (i in 1:nrHMM) {
		ind = which(subtrees[,2] == i)
		read.names = subtrees[ind, 1]
		avail.names = names(transcript.fas)
		subtree[[i]] = transcript.fas[intersect(read.names, avail.names)]
		sum = sum + length(subtree[[i]])
		if(sum(is.na(names(subtree[[i]]))) != 0) stop(" Fetching subset of fasta file failed.\n\n")
		# save the fasta for the subtrees
		subtree.fn[[i]] = paste(project.id, "_subtree_", i, ".fasta", sep="")
		out.fas = file.path(getwd(), outfolder, subtree.fn[[i]])
		write.fasta(subtree[[i]], names = names(subtree[[i]]), file.out = out.fas, open = "w", nbchar = 60, as.string = TRUE)
	}
	if(sum != nr.seqs) {
		cat(" WARNING: SUBTREEing did not pick up all reads.\n\n")
		put.warning(outhtml, paste("SUBTREEing did not pick up all reads listed in '", subtrees.txt, "'"))
	}
	if(length(subtree) != nrHMM) stop(" Not all subtrees found.\n\n")
	if(length(subtree.fn) != nrHMM) stop(" Not all subtrees found.\n\n")

	for (i in 1:nrHMM) {
		insert.blanks(outhtml, 2)
		insert.text(outhtml, text=paste("  Subtree", i, ":"), color="black")
		fn = file.path(getwd(), outfolder, subtree.fn[[i]])
		set.link(outhtml, fn, caption="fasta")
		insert.blanks(outhtml,2)
		insert.text(outhtml, text=paste("(", return.stained(length(subtree[[i]]), "blue"), "sequences )"), color="black")
		line.feed(outhtml, 1)
	}

} else {

	cat(paste("  Subtreeing is conducted by clustering.\n\n"))
	# insert.text(outhtml, text="Subtreeing is conducted by clustering.", color="black")
	# line.feed(outhtml, 1)
    #
	# if(exists("phylip.all") && (phylip.all != "")) {
	# 	if(file.exists(phylip.all)) {
	# 		cat(paste("  Alignment for all sequences will be loaded from '", phylip.all, "' (specified in configuration file)\n"))
	# 		aln <- read.alignment(phylip.all, format="phylip", forceToLower = FALSE)
	# 		if(class(aln) != "alignment") stop(paste(" Could not read alignment '", out.phy, "'.\n\n"))
	# 		if(aln$nb != length(transcript.fas))
	# 			stop(paste(" Alignment has", aln$nb, "entries, but transcript had", length(transcript.fas), "entries.\n\t Did you specify the correct alignment file?\n\n"))
	# 		success <- file.copy(phylip.all, file.path(getwd(), outfolder), overwrite = TRUE)
	# 		if(!success) stop(paste(" File '", phylip.all, "' could not be copied to folder", outfolder, ".\n\n"))
	# 		insert.text(outhtml, text="Loaded alignment for transcript sequence:", color="black")
	# 		set.link(outhtml, file.path(getwd(), outfolder, phylip.all), caption="phylip format")
	# 		line.feed(outhtml, 1)
	# 	} else {
	# 		stop(paste(" Alignment file '", phylip.all, "' was specified in configuration file, but was not found in the current folder.\n\n"))
	# 	}
	# } else {
	# 	fn = paste(project.id, "align1.phy", sep="_")
	# 	out.phy = file.path(getwd(), outfolder, fn)
	# 	clustal.log = paste(project.id, "align1_clustal.log", sep="_")
	# 	fn = file.path(getwd(), outfolder, clustal.log)
	# 	aln.res <- align.clustal(start.fasta.location, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = fn)
	# 	if(!file.exists(out.phy))  stop(" Could not complete clustalw alignment.\n\n")
	# 	aln <- read.alignment(out.phy, format="phylip", forceToLower = FALSE)
	# 	if(class(aln) != "alignment") stop(paste(" Could not read alignment '", out.phy, "'.\n\n"))
	# 	if(aln$nb != length(transcript.fas)) stop(paste(" Alignment has", aln$nb, "entries, but transcript had", length(transcript.fas), "entries.\n\n"))
	# 	insert.text(outhtml, text="Calculated alignment for transcript sequence:", color="black")
	# 	set.link(outhtml, out.phy, caption="phylip format")
	# 	line.feed(outhtml, 1)
	# }
    #
	# aln.1 = as.phyDat(aln)
	# if(class(aln.1) != "phyDat") stop(" Could not convert alignment to phyDat format.\n\n")
    #
	# distmat = dist.ml(aln.1, model="JC69")
	# if(class(distmat) != "dist") stop(" Could not calculate distance matrix from alignment.\n\n")
	# if(sum(is.na(distmat))	!= 0) stop(" Distance matrix includes missing data.\n\n")
    #
	# cluster = hclust(distmat, method = "ward.D")
	# if(class(cluster) != "hclust") stop(" Hierarchical clustering failed.\n\n")
    #
	# if(length(cluster$labels) < 25) {
	# 	labels = gsub(" ", "", cluster$labels)
	# } else {
	# 	labels = FALSE
	# }
	# if(is.logical(labels)) cat("  Labels in the phylogenetic trees cannot be shown: too many sequences.\n\n")
    #
	# X11(width=6, height=6)
	# plot(cluster, labels = labels, hang = 0.1, ann = TRUE, main = paste("Phylogenetic tree, project", project.id), sub = "", xlab = "", ylab = "", font.main=1)
	# mtext(side=1, paste(length(transcript.fas), "sequences"), col="blue", cex=0.9)
	# # fn = paste(project.id, "hc_tree_1.png", sep="_")
	# # outpng3 = file.path(getwd(), outfolder, fn)
	# # invisible(dev.copy(png, file = outpng3, unit="px", width=960, height=960, res=120)); invisible(dev.off())
	# fn = paste(project.id, "hc_tree_1.pdf", sep="_")
	# outpdf3 = file.path(getwd(), outfolder, fn)
	# invisible(dev.copy2pdf(file = outpdf3, out.type = "pdf")); invisible(dev.off())
	# cat(paste("  Phylogenetic tree (version 1) saved to", fn, "\n"))
	# # link plot
	# insert.text(outhtml, text="Phylogenetic tree, rectangle format (hclust):", color="black")
	# set.link(outhtml, outpdf3, caption="pdf")
	# insert.blanks(outhtml, 2)
	# # link tree in Newick format:
	# cluster.phylo <- as.phylo(cluster)
	# fn = "tree1.new"
	# tree.file = file.path(getwd(), outfolder, fn)
    # write.tree(phy = cluster.phylo, file = tree.file)
    # cat(paste("  Neighbor joining tree (version 1) saved to '", fn, "'\n\n"))
	# set.link(outhtml, tree.file, caption="Newick")
	# line.feed(outhtml, 1)
    #
    # cluster2 = nj(distmat)
	# if(class(cluster2) != "phylo") stop(" Could not calculate phylogenetic tree by neighbor-joining.\n\n")
    #
    # if(cluster2$Nnode < 25) {
	# 	labels = TRUE
	# } else {
	# 	labels = FALSE
	# }
    #
	# X11(width=6, height=6)
	# plot.phylo(cluster2, type = "unrooted", show.tip.label = labels, main = paste("Phylogenetic tree, project", project.id), font.main = 1)
	# fn = paste(project.id, "hc_tree_2.png", sep="_")
	# outpdf7 = file.path(getwd(), outfolder, fn)
	# invisible(dev.copy2pdf(file = outpdf7, out.type = "pdf")); invisible(dev.off())
	# cat(paste("  Phylogenetic tree (version 2) saved to", fn, "\n"))
	# # link plot
	# insert.text(outhtml, text="Phylogenetic tree, radial format (nj, ape):", color="black")
	# set.link(outhtml, outpdf7, caption="pdf")
	# insert.blanks(outhtml, 2)
	# # link tree in Newick format:
	# fn = "tree2.new"
	# tree.file = file.path(getwd(), outfolder, fn)
    # write.tree(phy = cluster2, file = tree.file)
    # cat(paste("  Neighbor joining tree (version 2) saved to '", fn, "'\n\n"))
    # set.link(outhtml, tree.file, caption="Newick")
	# line.feed(outhtml, 1)
    #
	# options(warn=-1)
	# repeat {
	# 	if(interactive())  {
	# 		nrHMM <- readline(prompt="  Number of Hidden Markov Models to create?(integer): ")
	# 	} else {
	# 		cat("  Number of Hidden Markov Models to create?(integer): ")
	# 		nrHMM = readLines("stdin", n=1)
	# 	}
	# 	nrHMM = as.integer(nrHMM)
	# 	if(!is.na(nrHMM))  {
	# 		if(nrHMM == round(nrHMM)) break
	# 	}
	# }
	# options(warn=0)
    #
	# cat(paste("\n  Building", nrHMM, "profile Hidden Markov Models:\n"))
	# insert.line(outhtml, percent=100, size=1)
	# insert.text(outhtml, text=paste("Number of HMM models to build:", return.stained(nrHMM, "red")), color="black")
	# line.feed(outhtml, 1)

	sub.trees = cutree(cluster, k = nrHMM)
	if(class(sub.trees) != "integer") stop(" Cutting into subtrees did not work.\n\n")
	new.names = gsub(" ", "", names(sub.trees))
	names(sub.trees) = new.names

	subtree = list()
	subtree.fn = list()
	sum = 0
	for (i in 1:nrHMM) {
		ind = which(sub.trees == i)
		read.names = names(sub.trees[ind])
		avail.names = names(transcript.fas)
		subtree[[i]] = transcript.fas[intersect(read.names, avail.names)]
		sum = sum + length(subtree[[i]])
		if(sum(is.na(names(subtree[[i]]))) != 0) stop(" Fetching subset of fasta file failed.\n\n")
		subtree.fn[[i]] = paste(project.id, "_subtree_", i, ".fasta", sep="")
		out.fas = file.path(getwd(), outfolder, subtree.fn[[i]])
		write.fasta(subtree[[i]], names = names(subtree[[i]]), file.out = out.fas, open = "w", nbchar = 60, as.string = TRUE)
	}
	if(sum != nr.seqs) {
		cat(" WARNING: SUBTREEing did not pick up all reads.\n\n")
		put.warning(outhtml, paste("SUBTREEing did not pick up all reads listed in '", subtrees.txt, "'"))
	}

	if(length(subtree) != nrHMM) stop(" Not all subtrees found (subtreeing by clustering).\n\n")
	if(length(subtree.fn) != nrHMM) stop(" Not all subtrees found (subtreeing by clustering).\n\n")

	# # link the results on website:
	# for (i in 1:nrHMM) {
	# 	insert.blanks(outhtml, 2)
	# 	insert.text(outhtml, text=paste("  Subtree", i, ":"), color="black")
	# 	fn = file.path(getwd(), outfolder, subtree.fn[[i]])
	# 	set.link(outhtml, fn, caption="fasta")
	# 	insert.blanks(outhtml, 2)
	# 	insert.text(outhtml, text=paste("(", length(subtree[[i]]), "sequences" , ")"), color="black")
	# 	line.feed(outhtml, 1)
	}

}

sum = 0
for (i in 1:nrHMM) sum = sum + length(subtree[[i]])
if(sum != length(transcript.fas)) stop(" Subtreeing did not collect all sequences).\n\n")

sizes = unlist(lapply(subtree, length))
if(min(sizes) < 3) {
	ind.to.remove = which(sizes < 3)
	for (i in ind.to.remove) subtree[[i]] <- NULL
	for (i in ind.to.remove) subtree.fn[[i]] <- NULL
	cat(paste("\n  The following subtrees have been removed because they contain less than 3 sequences:", paste(" Subtree", ind.to.remove, collapse=" , "), "\n\n"))
	nrHMM = length(subtree)
	# insert.line(outhtml, percent=100, size=1)
	warn.txt = paste("  The following subtrees have been removed because they contain less than 2 sequences:", paste(" Subtree", ind.to.remove, collapse=" , "))
	# put.warning(outhtml, warn.txt)
	if(nrHMM == 0) stop("  No subtrees left. Exiting.\n\n")
}



## +++ Create multiple alignment for each of the subtrees:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text=paste("Calculation of the alignments for the", nrHMM, "subtrees:"), color="black")
# line.feed(outhtml, 1)

subtree.alignment = list()
subtree.alignment.fn = list()
for (i in 1:nrHMM) {
	subtree.alignment.fn[[i]] = paste(project.id, "_subtree_", i, ".phy", sep="")
	out.phy = file.path(getwd(), outfolder, subtree.alignment.fn[[i]])
	in.fas = file.path(getwd(), outfolder, subtree.fn[[i]])
	clustal.log = paste(project.id, "_subtree_", i, "_clustalw.log", sep="")
	fn = file.path(getwd(), outfolder, clustal.log)
	aln.res <- align.clustal(in.fas, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = fn, nr.cpu = clustalo.nr.cpu)
	if(!file.exists(out.phy))  stop(paste(" Could not complete clustalw alignment (step2) for subtree", i, ".\n\n"))
	subtree.alignment[[i]] <- read.alignment(out.phy, format="phylip", forceToLower = FALSE)
	if(class(subtree.alignment[[i]]) != "alignment") stop(paste(" Could not read alignment '", subtree.alignment.fn[[i]], "' (clustalw step 2).\n\n"))
	if(subtree.alignment[[i]]$nb != length(subtree[[i]])) stop(paste(" Alignment has", subtree.alignment[[i]]$nb, "entries, but transcript had", length(transcript.fas), "entries.\n\n"))
	# insert.blanks(outhtml, 2)
	# insert.text(outhtml, text=paste("Calculated alignment for subtree", i, ":"), color="black")
	# set.link(outhtml, out.phy, caption="phylip format")
	# line.feed(outhtml, 1)
}
cat("\n")
if(length(subtree.alignment) != nrHMM) stop(" Not all subtrees found (aligning).\n\n")
if(length(subtree.alignment.fn) != nrHMM) stop(" Not all subtrees found (aligning).\n\n")


## +++ Calculate distance matrix for each of the subtree alignments:

subtree.alignment.1 = list()
subtree.distance = list()
for (i in 1:nrHMM) {
	subtree.alignment.1[[i]] = as.phyDat(subtree.alignment[[i]])
	if(class(subtree.alignment.1[[i]])	!= "phyDat") stop(" Conversion to phyDat format not successful.\n\n")
	subtree.distance[[i]] = dist.ml(subtree.alignment.1[[i]], model="JC69")
    if(class(subtree.distance[[i]])	!= "dist")  stop(" Subtree distance calculation not successful.\n\n")
	if(sum(is.na(subtree.distance[[i]])) != 0) stop(" Subtree distance calculation produced missing data.\n\n")
}
if(length(subtree.alignment.1) != nrHMM) stop(" Not all subtree distances calculated.\n\n")
if(length(subtree.distance) != nrHMM) stop(" Not all subtree distances calculated.\n\n")


## +++ Conduct hierarchical clustering:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Dendograms for the subtrees:", color="black")
# line.feed(outhtml, 1)

subtree.cluster = list()
for (i in 1:nrHMM) {
	subtree.cluster[[i]] = hclust(subtree.distance[[i]], method = "complete")
	if(class(subtree.cluster[[i]]) != "hclust") stop(" Clustering of subtree not successful.\n\n")
   	if(length(subtree.cluster[[i]]$labels) < 25) {
		labels = gsub(" ", "", subtree.cluster[[i]]$labels)
	} else {
		labels = FALSE
	}
	if(is.logical(labels)) cat("  Labels in the subtree cannot be shown: too many sequences.\n")

	fn = paste("subtree_cluster_", i, ".png", sep="")
	fn = file.path(getwd(), outfolder, fn)
    png(filename = fn, units = "px", width = 960, height = 960, res=120)
    maintxt = paste("Phylogenetic tree, project", project.id, ", subtree", i)
	plot(subtree.cluster[[i]], labels = labels, hang = 0.1, ann = TRUE, main = maintxt, sub = "", xlab = "", ylab = "", font.main=1)
    dev.off()
	# insert.blanks(outhtml, 2)
	# insert.text(outhtml, text=paste("Phylogenetic tree for subtree", i, ": "), color="black")
	# set.link(outhtml, fn, caption="png")
	# line.feed(outhtml, 1)
}


## +++ Cut subtree into appropiate number of subgroups:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Picking from subtrees ('Thinning'):", color="black")
# line.feed(outhtml, 1)

thinned = logical(length = nrHMM)
subtree.group = list()

picked.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
rownames(picked.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(picked.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

for (i in 1:nrHMM) {
    subtree.size = length(subtree.cluster[[i]]$labels)
    subtree.size2 = length(subtree[[i]])
    if(subtree.size != subtree.size2) stop(paste(" Fasta for subtree does not correspond to referring cluster: i =",i, "\n\n"))
	num.clust = max(round(subgroup.percent / 100 * subtree.size),1)

	if(num.clust == 1) {
		cat(paste("  Cluster", i, " - not subdividing - too small, only", subtree.size, "members\n"))
		fasta.picked = subtree[[i]]
		picked.fn[i,1] = paste("subtree", i, "picked.fasta", sep="_")
		fn = file.path(getwd(), outfolder, picked.fn[i,1])
		write.fasta(fasta.picked, names = names(fasta.picked), file.out = fn, open = "w", nbchar = 60, as.string = TRUE)
		# link on webpage:
		# insert.blanks(outhtml,4)
		# insert.text(outhtml, text=paste("Picked from subtree", i, ":"), color="black")
		# set.link(outhtml, fn, caption="fas")
		# insert.blanks(outhtml,2)
		# insert.text(outhtml, text=paste("( no trials made - too few sequences )"), color="black")
		# line.feed(outhtml, 1)
	} else {
		thinned[i] = TRUE
		cat(paste("  Cluster", i, " - subdividing into", num.clust, "subgroups\n"))
		subtree.group[[i]] = cutree(subtree.cluster[[i]], k = num.clust)
		if(class(subtree.group[[i]]) != "integer") stop(" Cutting of subtree not successful.\n\n")
		if(length(subtree[[i]]) != length(subtree.cluster[[i]]$labels)) stop(paste(" subtree object does not refer to subtree.cluster object: i =",i, "\n\n"))
		if(length(subtree[[i]]) != length(subtree.group[[i]])) stop(paste(" subtree object does not refer to subtree.group object: i =",i, "\n\n"))

		fn = paste("subtree_cluster_", i, "_groupsize.png", sep="")
		fn = file.path(getwd(), outfolder, fn)
		png(filename = fn, units = "px", width = 960, height = 960, res=120)
		barplot(table(subtree.group[[i]]), col="red", main=paste("Subgroup sizes for subtree",i), font.main=1)
		mtext(side=3, paste(project.id, "subtree", i, ",", subtree.size, "sequences,", num.clust, "subgroups"), col="blue", cex=0.9)
		dev.off()
		# insert.blanks(outhtml, 2)
		# insert.text(outhtml, text=paste("Number of sequences in the groups of subtree", i, ": "), color="black")
		# set.link(outhtml, fn, caption="png")
		# line.feed(outhtml, 1)
		to_pick = matrix(NA, nrow = num.clust, ncol = nr.trials.random.picking)
		# insert.blanks(outhtml,4)
		# insert.text(outhtml, text=paste("Picked from subtree", i, ":"), color="black")
		for (k in 1:nr.trials.random.picking) {
			to_pick[,k] = pick_one_from_each_cluster(subtree.group[[i]], silent = TRUE)
			fasta.picked = subtree[[i]][to_pick[,k]]
			if(length(fasta.picked) != num.clust) stop(paste(" Fasta incorrectly picked: i =",i, ", k =", k,"\n\n"))
			picked.fn[i,k] = paste("subtree", i, "trial", k, "picked.fasta", sep="_")
			fn = file.path(getwd(), outfolder, picked.fn[i,k])
			write.fasta(fasta.picked, names = names(fasta.picked), file.out = fn, open = "w", nbchar = 60, as.string = TRUE)
			# link on webpage:
			# insert.blanks(outhtml,1)
			# insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
			# set.link(outhtml, fn, caption="fas")
			# if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
		}
		colnames(to_pick) = paste("trial", 1:nr.trials.random.picking, sep="_")
		rownames(to_pick) = paste("cl", 1:num.clust, sep="_")
		if(sum(is.na(to_pick)) != 0) stop("  Could not find indices for sequences to pick from subgroups.\n\n")
	}
}
cat("\n")


## Create profile HMM from each of the fasta files:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Create profile Hidden Markov Models for subtrees:", color="black")
# line.feed(outhtml, 1)

phyfile.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
rownames(phyfile.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(phyfile.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

for (i in 1:nrHMM) {
	if(!thinned[i]) {
		fas.fn = picked.fn[i,1]
		fas.fn = file.path(getwd(), outfolder, fas.fn)
		if(!file.exists(fas.fn)) stop("  Fasta file not found for alignment: i = ", i, " (no thinning for this subgroup, no trials).\n\n")
		phyfile.fn[i, 1] = paste("subtree", i, "no_trial_aln.phy", sep="_")
		out.phy = file.path(getwd(), outfolder, phyfile.fn[i, 1])
		clustal.log = paste("subtree", i, "no_trial_clustal.log", sep="_")
		log.fn = file.path(getwd(), outfolder, clustal.log)
		aln.res <- align.clustal(fas.fn, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = log.fn)
		if(!file.exists(out.phy))  stop(paste(" Could not complete clustalw alignment for subtree", i, " (no thinning for this subgroup, no trials).\n\n"))
	} else {
		for (k in 1:nr.trials.random.picking) {
			fas.fn = picked.fn[i,k]
			fas.fn = file.path(getwd(), outfolder, fas.fn)
			if(!file.exists(fas.fn)) stop("  Fasta file not found for alignment: i = ", i, ",  k = ", k, "\n\n")
			# define output file (phylip format):
			phyfile.fn[i, k] = paste("subtree", i, "trial", k, "aln.phy", sep="_")
			out.phy = file.path(getwd(), outfolder, phyfile.fn[i, k])
			clustal.log = paste("subtree", i, "trial", k, "clustal.log", sep="_")
			log.fn = file.path(getwd(), outfolder, clustal.log)
			aln.res <- align.clustal(fas.fn, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = log.fn)
			if(!file.exists(out.phy))  stop(paste(" Could not complete clustalw alignment for subtree", i, ", trial", k, ".\n\n"))
		}
	}
}

hmmfile.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
rownames(hmmfile.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(hmmfile.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

cat("\n")
for (i in 1:nrHMM) {
	if(!thinned[i]) {
		phy.fn = phyfile.fn[i, 1]
		phy.fn = file.path(getwd(), outfolder, phy.fn)
		if(!file.exists(phy.fn)) stop(paste("  Phylip file '", phy.fn,  "' not found.\n\n"))
		hmmfile.fn[i, 1] = paste("subtree", i, "no_trial.hmm", sep="_")
		out.hmm = file.path(getwd(), outfolder, hmmfile.fn[i, 1])
		HMM.name = paste("HMM_subtree", i, "no_trial", sep="_")
		summary.file = paste(HMM.name, "summary.txt", sep = "_")
		summary.file = file.path(getwd(), outfolder, summary.file)
		new.alignment = paste(HMM.name, "new_align.txt", sep = "_")
		new.alignment = file.path(getwd(), outfolder, new.alignment)
		hmmbuild.res <- build.hmm(alignment = phy.fn, outhmm = out.hmm, hmm.name = HMM.name, summary.file = summary.file, new.alignment = new.alignment, informat = "PHYLIP", nr.cpu = hmmbuild.nr.cpu)
		if(!file.exists(out.hmm))  stop(paste(" Could not complete hmmbuild for subtree", i, " (no thinning for this subgroup, no trials).\n\n"))
		# link on webpage:
		# insert.blanks(outhtml,2)
		# insert.text(outhtml, text=paste("HMM for subtree", i, ":"), color="black")
		# set.link(outhtml, out.hmm, caption="txt")
		# insert.blanks(outhtml,2)
		# insert.text(outhtml, text=paste("( no trials made - too few sequences )"), color="black")
		# line.feed(outhtml, 1)
	} else {
		# insert.blanks(outhtml,2)
		# insert.text(outhtml, text=paste("HMM's for subtree", i, ":"), color="black")
		for (k in 1:nr.trials.random.picking) {
			phy.fn = phyfile.fn[i, k]
			phy.fn = file.path(getwd(), outfolder, phy.fn)
			if(!file.exists(phy.fn)) stop(paste("  Phylip file '", phy.fn,  "' not found.\n\n"))
			hmmfile.fn[i, k] = paste("subtree", i, "trial", k, "profile.hmm", sep="_")
			out.hmm = file.path(getwd(), outfolder, hmmfile.fn[i, k])
			HMM.name = paste("HMM_subtree", i, "trial", k, sep="_")
			summary.file = paste(HMM.name, "summary.txt", sep = "_")
			summary.file = file.path(getwd(), outfolder, summary.file)
			new.alignment = paste(HMM.name, "new_align.txt", sep = "_")
			new.alignment = file.path(getwd(), outfolder, new.alignment)
			hmmbuild.res <- build.hmm(alignment = phy.fn, outhmm = out.hmm, hmm.name = HMM.name, summary.file = summary.file, new.alignment = new.alignment, informat = "PHYLIP", nr.cpu = hmmbuild.nr.cpu)
			if(!file.exists(out.hmm))  stop(paste(" Could not complete ' hmmbuild'  for subtree", i, ", trial", k, ".\n\n"))
			# link on webpage:
			# insert.blanks(outhtml,1)
			# insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
			# set.link(outhtml, out.hmm, caption="txt")
			# if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
		}
	}
}


## +++ zip HMMer output files
param.file = "quantizyme_model_params.RData"
params = list()
params[["nrHMM"]] = nrHMM
params[["thinned"]] = thinned
params[["hmmfile.fn"]] = hmmfile.fn
params[["nr.trials.random.picking"]] = nr.trials.random.picking
save(params, file = param.file)
cres = file.copy(param.file, outfolder)
if(!cres) stop(paste(" Could not copy file '", param.file, "' to", outfolder, "\n\n"))
current.location = getwd()
setwd(paste(current.location,  outfolder, sep = "/"))
zipfile = paste(project.id, "_models.zip", sep = "")
models.names = as.vector(hmmfile.fn)
models.names = models.names[!is.na(models.names)]
number.models = length(models.names)
models.names = paste(models.names, collapse = " ")
models.names = paste(models.names, param.file, collapse = " ")
command = paste("zip -q", zipfile,  models.names)
cat(paste("\n  Zipping models and parameters\n"))
zipres <- try(system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL))
if(zipres != 0) stop("\n\n  Zipping the Hidden Markov Models failed.\n\n")
zip.list = paste(project.id, "HMM_listing.txt", sep = "_")
c2 = paste("unzip -l", zipfile, ">", zip.list)
listres <- try(system(c2, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL))
if(listres != 0) stop("\n\n  Listing the Hidden Markov Models failed.\n\n")
setwd(current.location)
# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text = "Hidden Markov Models:", color = "black")
zip.location = file.path(getwd(), outfolder, zipfile)
# line.feed(outhtml, 1)
# insert.blanks(outhtml,2)
# insert.text(outhtml, text="Models:", color="black")
# insert.blanks(outhtml,1)
# set.link(outhtml, zip.location, caption = zipfile)
# insert.text(outhtml, text=paste("(", return.stained(number.models, "blue"), "models )"), color="black")
# line.feed(outhtml, 1)
# insert.blanks(outhtml,2)
# insert.text(outhtml, text="Listing of models + parameter file:", color="black")
# insert.blanks(outhtml,1)
list.location = file.path(getwd(), outfolder, zip.list)
# set.link(outhtml, list.location, caption = zip.list)
# insert.blanks(outhtml, 1)
# line.feed(outhtml, 1)


## +++ End
# insert.line(outhtml, percent=100, size=1)
# stop.html(outhtml)
# cat(paste("\n  Results written to '", outhtml, "'\n\n"))

# uwe.menzel@matstat.de
# uwe.menzel@gmail.com
