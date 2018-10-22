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
    make_option(c("-d", "--outdir"), type="character", default=NULL,
            help="outdir [default= %default] (MANDATORY)"),
    make_option(c("-e", "--environment_file"), type="character", default=NULL,
        help="R environment file [default= %default] (MANDATORY)", metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="input fasta file [default= %default] (MANDATORY)"),
    make_option(c("-n", "--n_subtrees"), type="integer", default=NULL,
        help="number of subtrees [default= %default] (MANDATORY)", metavar="character"),
    make_option(c("-p", "--projectID"), type="character", default=NULL,
            help="projectID [default= %default] (MANDATORY)"),
    make_option(c("-s", "--subtree_percent"), type="integer", default=NULL,
            help="subtree_percent_group [default= %default] (MANDATORY)")
    make_option(c("-z", "--done_file"), type="integer", default=NULL,
            help="flag file [default= %default] (MANDATORY)")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

out.phy <- opt$a
envcluster_file = opt$e
nrHMM = as.integer(opt$n)
transcript.fasta = opt$i
project.id = opt$p
outfolder = opt$outdir
i = opt$n
subgroup.percent = opt$s
flag_file = opt$z

pick_one_from_each_cluster <- function(clusters, silent = FALSE) {
	## uwe.menzel@slu.se  uwe.menzel@gmail.com
	if(class(clusters) != "integer") stop("  pick_one_from_each_cluster: Input must be of type 'integer'.\n\n")
	num.cl = length(unique(clusters))
	picked = integer(num.clust)

	for (i in 1:num.cl) {
		if(!silent) cat(paste("Cluster", i, "\n"))
		indx = which(clusters == i)
		if(length(indx) == 1) {
			if(!silent) cat(paste("  Available:", as.numeric(indx), "\n"))
			ix = indx
			if(!silent) cat(paste("  Choosen:  ", as.numeric(ix), "\n\n"))
		} else {
			if(!silent) cat(paste("  Available:", paste(as.numeric(indx), collapse=" "), "\n", collapse=" "))
			ix = sample(indx, 1)
			if(!silent) cat(paste("  Choosen:  ", as.numeric(ix), "\n\n"))
		}
		picked[i] = ix
	}
	return(picked)
}

## +++ Calculate distance matrix for each of the subtree alignments:

aln <- read.alignment(out.phy, format="phylip", forceToLower = FALSE)
if(class(aln) != "alignment") stop(paste(" Could not read alignment '", out.phy, "'.\n\n"))
if(aln$nb != length(transcript.fas)) stop(paste(" Alignment has", aln$nb, "entries, but transcript had", length(transcript.fas), "entries.\n\n"))

aln.1 = as.phyDat(aln)
if(class(aln.1)	!= "phyDat") stop(" Conversion to phyDat format not successful.\n\n")
aln.1.distance = dist.ml(aln.1, model="JC69")
if(class(aln.1.distance)	!= "dist")  stop(" Subtree distance calculation not successful.\n\n")
if(sum(is.na(aln.1.distance)) != 0) stop(" Subtree distance calculation produced missing data.\n\n")

# subtree.alignment.1 = list()
# subtree.distance = list()
# for (i in 1:nrHMM) {
# 	subtree.alignment.1[[i]] = as.phyDat(subtree.alignment[[i]])
# 	if(class(subtree.alignment.1[[i]])	!= "phyDat") stop(" Conversion to phyDat format not successful.\n\n")
# 	subtree.distance[[i]] = dist.ml(subtree.alignment.1[[i]], model="JC69")
#     if(class(subtree.distance[[i]])	!= "dist")  stop(" Subtree distance calculation not successful.\n\n")
# 	if(sum(is.na(subtree.distance[[i]])) != 0) stop(" Subtree distance calculation produced missing data.\n\n")
# }
# if(length(subtree.alignment.1) != nrHMM) stop(" Not all subtree distances calculated.\n\n")
# if(length(subtree.distance) != nrHMM) stop(" Not all subtree distances calculated.\n\n")


## +++ Conduct hierarchical clustering:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Dendograms for the subtrees:", color="black")
# line.feed(outhtml, 1)

# subtree.cluster = list()
# for (i in 1:nrHMM) {
# 	subtree.cluster[[i]] = hclust(subtree.distance[[i]], method = "complete")
# 	if(class(subtree.cluster[[i]]) != "hclust") stop(" Clustering of subtree not successful.\n\n")
#    	if(length(subtree.cluster[[i]]$labels) < 25) {
# 		labels = gsub(" ", "", subtree.cluster[[i]]$labels)
# 	} else {
# 		labels = FALSE
# 	}
# 	if(is.logical(labels)) cat("  Labels in the subtree cannot be shown: too many sequences.\n")
#
# 	fn = paste("subtree_cluster_", i, ".png", sep="")
# 	fn = file.path(getwd(), outfolder, fn)
#     png(filename = fn, units = "px", width = 960, height = 960, res=120)
#     maintxt = paste("Phylogenetic tree, project", project.id, ", subtree", i)
# 	plot(subtree.cluster[[i]], labels = labels, hang = 0.1, ann = TRUE, main = maintxt, sub = "", xlab = "", ylab = "", font.main=1)
#     dev.off()
# 	# insert.blanks(outhtml, 2)
# 	# insert.text(outhtml, text=paste("Phylogenetic tree for subtree", i, ": "), color="black")
# 	# set.link(outhtml, fn, caption="png")
# 	# line.feed(outhtml, 1)
# }

subtree.cluster = hclust(aln.1.distance, method = "complete")
if(class(subtree.cluster) != "hclust") stop(" Clustering of subtree not successful.\n\n")
if(length(subtree.cluster$labels) < 25) {
    labels = gsub(" ", "", subtree.cluster$labels)
} else {
    labels = FALSE
}
if(is.logical(labels)) cat("  Labels in the subtree cannot be shown: too many sequences.\n")

fn = paste("subtree_cluster_", i, ".png", sep="")
fn = file.path(outfolder, fn)
png(filename = fn, units = "px", width = 960, height = 960, res=120)
maintxt = paste("Phylogenetic tree, project", project.id, ", subtree", i)
plot(subtree.cluster, labels = labels, hang = 0.1, ann = TRUE, main = maintxt, sub = "", xlab = "", ylab = "", font.main=1)
dev.off()
# insert.blanks(outhtml, 2)
# insert.text(outhtml, text=paste("Phylogenetic tree for subtree", i, ": "), color="black")
# set.link(outhtml, fn, caption="png")
# line.feed(outhtml, 1)


## +++ Cut subtree into appropiate number of subgroups:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Picking from subtrees ('Thinning'):", color="black")
# line.feed(outhtml, 1)

thinned = logical(length = nrHMM)
subtree.group = list()

picked.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
rownames(picked.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(picked.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

subtree.size = length(subtree.cluster$labels)
# subtree.size2 = length(subtree[[i]])
# if(subtree.size != subtree.size2) stop(paste(" Fasta for subtree does not correspond to referring cluster: i =",i, "\n\n"))
num.clust = max(round(subgroup.percent / 100 * subtree.size),1)

phyfile.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
rownames(phyfile.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(phyfile.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

hmmfile.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
rownames(hmmfile.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(hmmfile.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

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
    fas.fn = picked.fn[i,1]
    fas.fn = file.path(getwd(), outfolder, fas.fn)
    if(!file.exists(fas.fn)) stop("  Fasta file not found for alignment: i = ", i, " (no thinning for this subgroup, no trials).\n\n")
    phyfile.fn[i, 1] = paste("subtree", i, "no_trial_aln.phy", sep="_")
    out.phy = file.path(getwd(), outfolder, phyfile.fn[i, 1])
    clustal.log = paste("subtree", i, "no_trial_clustal.log", sep="_")
    log.fn = file.path(getwd(), outfolder, clustal.log)
    aln.res <- align.clustal(fas.fn, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = log.fn)
    if(!file.exists(out.phy))  stop(paste(" Could not complete clustalw alignment for subtree", i, " (no thinning for this subgroup, no trials).\n\n"))

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
	colnames(to_pick) = paste("trial", 1:nr.trials.random.picking, sep="_")
	rownames(to_pick) = paste("cl", 1:num.clust, sep="_")
	if(sum(is.na(to_pick)) != 0) stop("  Could not find indices for sequences to pick from subgroups.\n\n")
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
		# insert.blanks(outhtml,2)
		# insert.text(outhtml, text=paste("HMM's for subtree", i, ":"), color="black")
	for (k in 1:nr.trials.random.picking) {
		phy.fn = phyfile.fn[i, k]
		phy.fn = file.path(getwd(), outfolder, phy.fn)
		if(!file.exists(phy.fn)) stop(paste("  Phylip file '", phy.fn,  "' not found.\n\n"))
		hmmfile.fn[i, k] = paste("subtree", i, "trial", k, "profile.hmm", sep="_")
		out.hmm = file.path(outfolder, hmmfile.fn[i, k])
		HMM.name = paste("HMM_subtree", i, "trial", k, sep="_")
		summary.file = paste(HMM.name, "summary.txt", sep = "_")
		summary.file = file.path(outfolder, summary.file)
		new.alignment = paste(HMM.name, "new_align.txt", sep = "_")
		new.alignment = file.path(outfolder, new.alignment)
		hmmbuild.res <- build.hmm(alignment = phy.fn, outhmm = out.hmm, hmm.name = HMM.name, summary.file = summary.file, new.alignment = new.alignment, informat = "PHYLIP", nr.cpu = hmmbuild.nr.cpu)
		if(!file.exists(out.hmm))  stop(paste(" Could not complete ' hmmbuild'  for subtree", i, ", trial", k, ".\n\n"))
		# link on webpage:
		# insert.blanks(outhtml,1)
		# insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
		# set.link(outhtml, out.hmm, caption="txt")
		# if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
	}
	# colnames(to_pick) = paste("trial", 1:nr.trials.random.picking, sep="_")
	# rownames(to_pick) = paste("cl", 1:num.clust, sep="_")
	# if(sum(is.na(to_pick)) != 0) stop("  Could not find indices for sequences to pick from subgroups.\n\n")
}
cat("\n")

file.create(flag_file, overwrite=TRUE)

#######
# for (i in 1:nrHMM) {
#     subtree.size = length(subtree.cluster[[i]]$labels)
#     subtree.size2 = length(subtree[[i]])
#     if(subtree.size != subtree.size2) stop(paste(" Fasta for subtree does not correspond to referring cluster: i =",i, "\n\n"))
# 	num.clust = max(round(subgroup.percent / 100 * subtree.size),1)
#
# 	if(num.clust == 1) {
# 		cat(paste("  Cluster", i, " - not subdividing - too small, only", subtree.size, "members\n"))
# 		fasta.picked = subtree[[i]]
# 		picked.fn[i,1] = paste("subtree", i, "picked.fasta", sep="_")
# 		fn = file.path(getwd(), outfolder, picked.fn[i,1])
# 		write.fasta(fasta.picked, names = names(fasta.picked), file.out = fn, open = "w", nbchar = 60, as.string = TRUE)
# 		# link on webpage:
# 		# insert.blanks(outhtml,4)
# 		# insert.text(outhtml, text=paste("Picked from subtree", i, ":"), color="black")
# 		# set.link(outhtml, fn, caption="fas")
# 		# insert.blanks(outhtml,2)
# 		# insert.text(outhtml, text=paste("( no trials made - too few sequences )"), color="black")
# 		# line.feed(outhtml, 1)
# 	} else {
# 		thinned[i] = TRUE
# 		cat(paste("  Cluster", i, " - subdividing into", num.clust, "subgroups\n"))
# 		subtree.group[[i]] = cutree(subtree.cluster[[i]], k = num.clust)
# 		if(class(subtree.group[[i]]) != "integer") stop(" Cutting of subtree not successful.\n\n")
# 		if(length(subtree[[i]]) != length(subtree.cluster[[i]]$labels)) stop(paste(" subtree object does not refer to subtree.cluster object: i =",i, "\n\n"))
# 		if(length(subtree[[i]]) != length(subtree.group[[i]])) stop(paste(" subtree object does not refer to subtree.group object: i =",i, "\n\n"))
#
# 		fn = paste("subtree_cluster_", i, "_groupsize.png", sep="")
# 		fn = file.path(getwd(), outfolder, fn)
# 		png(filename = fn, units = "px", width = 960, height = 960, res=120)
# 		barplot(table(subtree.group[[i]]), col="red", main=paste("Subgroup sizes for subtree",i), font.main=1)
# 		mtext(side=3, paste(project.id, "subtree", i, ",", subtree.size, "sequences,", num.clust, "subgroups"), col="blue", cex=0.9)
# 		dev.off()
# 		# insert.blanks(outhtml, 2)
# 		# insert.text(outhtml, text=paste("Number of sequences in the groups of subtree", i, ": "), color="black")
# 		# set.link(outhtml, fn, caption="png")
# 		# line.feed(outhtml, 1)
# 		to_pick = matrix(NA, nrow = num.clust, ncol = nr.trials.random.picking)
# 		# insert.blanks(outhtml,4)
# 		# insert.text(outhtml, text=paste("Picked from subtree", i, ":"), color="black")
# 		for (k in 1:nr.trials.random.picking) {
# 			to_pick[,k] = pick_one_from_each_cluster(subtree.group[[i]], silent = TRUE)
# 			fasta.picked = subtree[[i]][to_pick[,k]]
# 			if(length(fasta.picked) != num.clust) stop(paste(" Fasta incorrectly picked: i =",i, ", k =", k,"\n\n"))
# 			picked.fn[i,k] = paste("subtree", i, "trial", k, "picked.fasta", sep="_")
# 			fn = file.path(getwd(), outfolder, picked.fn[i,k])
# 			write.fasta(fasta.picked, names = names(fasta.picked), file.out = fn, open = "w", nbchar = 60, as.string = TRUE)
# 			# link on webpage:
# 			# insert.blanks(outhtml,1)
# 			# insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
# 			# set.link(outhtml, fn, caption="fas")
# 			# if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
# 		}
# 		colnames(to_pick) = paste("trial", 1:nr.trials.random.picking, sep="_")
# 		rownames(to_pick) = paste("cl", 1:num.clust, sep="_")
# 		if(sum(is.na(to_pick)) != 0) stop("  Could not find indices for sequences to pick from subgroups.\n\n")
# 	}
# }
# cat("\n")
########################

## Create profile HMM from each of the fasta files:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Create profile Hidden Markov Models for subtrees:", color="black")
# line.feed(outhtml, 1)

# phyfile.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
# rownames(phyfile.fn) = paste("subtree", 1:nrHMM, sep="_")
# colnames(phyfile.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

# for (i in 1:nrHMM) {
# 	if(!thinned[i]) {
# 		fas.fn = picked.fn[i,1]
# 		fas.fn = file.path(getwd(), outfolder, fas.fn)
# 		if(!file.exists(fas.fn)) stop("  Fasta file not found for alignment: i = ", i, " (no thinning for this subgroup, no trials).\n\n")
# 		phyfile.fn[i, 1] = paste("subtree", i, "no_trial_aln.phy", sep="_")
# 		out.phy = file.path(getwd(), outfolder, phyfile.fn[i, 1])
# 		clustal.log = paste("subtree", i, "no_trial_clustal.log", sep="_")
# 		log.fn = file.path(getwd(), outfolder, clustal.log)
# 		aln.res <- align.clustal(fas.fn, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = log.fn)
# 		if(!file.exists(out.phy))  stop(paste(" Could not complete clustalw alignment for subtree", i, " (no thinning for this subgroup, no trials).\n\n"))
# 	} else {
# 		for (k in 1:nr.trials.random.picking) {
# 			fas.fn = picked.fn[i,k]
# 			fas.fn = file.path(getwd(), outfolder, fas.fn)
# 			if(!file.exists(fas.fn)) stop("  Fasta file not found for alignment: i = ", i, ",  k = ", k, "\n\n")
# 			# define output file (phylip format):
# 			phyfile.fn[i, k] = paste("subtree", i, "trial", k, "aln.phy", sep="_")
# 			out.phy = file.path(getwd(), outfolder, phyfile.fn[i, k])
# 			clustal.log = paste("subtree", i, "trial", k, "clustal.log", sep="_")
# 			log.fn = file.path(getwd(), outfolder, clustal.log)
# 			aln.res <- align.clustal(fas.fn, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = log.fn)
# 			if(!file.exists(out.phy))  stop(paste(" Could not complete clustalw alignment for subtree", i, ", trial", k, ".\n\n"))
# 		}
# 	}
# }

# hmmfile.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)
# rownames(hmmfile.fn) = paste("subtree", 1:nrHMM, sep="_")
# colnames(hmmfile.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

cat("\n")
# for (i in 1:nrHMM) {
# 	if(!thinned[i]) {
# 		phy.fn = phyfile.fn[i, 1]
# 		phy.fn = file.path(getwd(), outfolder, phy.fn)
# 		if(!file.exists(phy.fn)) stop(paste("  Phylip file '", phy.fn,  "' not found.\n\n"))
# 		hmmfile.fn[i, 1] = paste("subtree", i, "no_trial.hmm", sep="_")
# 		out.hmm = file.path(getwd(), outfolder, hmmfile.fn[i, 1])
# 		HMM.name = paste("HMM_subtree", i, "no_trial", sep="_")
# 		summary.file = paste(HMM.name, "summary.txt", sep = "_")
# 		summary.file = file.path(getwd(), outfolder, summary.file)
# 		new.alignment = paste(HMM.name, "new_align.txt", sep = "_")
# 		new.alignment = file.path(getwd(), outfolder, new.alignment)
# 		hmmbuild.res <- build.hmm(alignment = phy.fn, outhmm = out.hmm, hmm.name = HMM.name, summary.file = summary.file, new.alignment = new.alignment, informat = "PHYLIP", nr.cpu = hmmbuild.nr.cpu)
# 		if(!file.exists(out.hmm))  stop(paste(" Could not complete hmmbuild for subtree", i, " (no thinning for this subgroup, no trials).\n\n"))
# 		# link on webpage:
# 		# insert.blanks(outhtml,2)
# 		# insert.text(outhtml, text=paste("HMM for subtree", i, ":"), color="black")
# 		# set.link(outhtml, out.hmm, caption="txt")
# 		# insert.blanks(outhtml,2)
# 		# insert.text(outhtml, text=paste("( no trials made - too few sequences )"), color="black")
# 		# line.feed(outhtml, 1)
# 	} else {
# 		# insert.blanks(outhtml,2)
# 		# insert.text(outhtml, text=paste("HMM's for subtree", i, ":"), color="black")
# 		for (k in 1:nr.trials.random.picking) {
# 			phy.fn = phyfile.fn[i, k]
# 			phy.fn = file.path(getwd(), outfolder, phy.fn)
# 			if(!file.exists(phy.fn)) stop(paste("  Phylip file '", phy.fn,  "' not found.\n\n"))
# 			hmmfile.fn[i, k] = paste("subtree", i, "trial", k, "profile.hmm", sep="_")
# 			out.hmm = file.path(getwd(), outfolder, hmmfile.fn[i, k])
# 			HMM.name = paste("HMM_subtree", i, "trial", k, sep="_")
# 			summary.file = paste(HMM.name, "summary.txt", sep = "_")
# 			summary.file = file.path(getwd(), outfolder, summary.file)
# 			new.alignment = paste(HMM.name, "new_align.txt", sep = "_")
# 			new.alignment = file.path(getwd(), outfolder, new.alignment)
# 			hmmbuild.res <- build.hmm(alignment = phy.fn, outhmm = out.hmm, hmm.name = HMM.name, summary.file = summary.file, new.alignment = new.alignment, informat = "PHYLIP", nr.cpu = hmmbuild.nr.cpu)
# 			if(!file.exists(out.hmm))  stop(paste(" Could not complete ' hmmbuild'  for subtree", i, ", trial", k, ".\n\n"))
# 			# link on webpage:
# 			# insert.blanks(outhtml,1)
# 			# insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
# 			# set.link(outhtml, out.hmm, caption="txt")
# 			# if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
# 		}
# 	}
# }
