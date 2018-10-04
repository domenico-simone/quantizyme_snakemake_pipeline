
# quantizyme_model.R    Draft 1.0   March 26, 2018  Uwe Menzel

# --------------------------------------------------------------------------------------------
#  R.Version()$platform			# "i686-pc-linux-gnu"
#  R.Version()$version.string	# "R version 3.4.2 (2017-09-28)"
#  getwd()						# "/home/uwe/Desktop/R/METAGENOMICS_PIPELINE"
#  date()						# "Mon Mar 26 12:41:43 2018"
# --------------------------------------------------------------------------------------------

# uwe.menzel@matstat.de
# uwe.menzel@gmail.com

## Rscript --vanilla  quantizyme_model.R   GT48  enzymes.fasta



## +++ Load R libraries and check if external programs are available:

required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs)
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))
suppressMessages(library(phangorn))
suppressMessages(library(ape))

system("echo $DISPLAY")
#Sys.setenv("DISPLAY"=":0.0")



# Check if required programs are installed:
clustalw.location <- Sys.which("clustalo")
if(clustalw.location == ""){
    print("Please install 'clustalo' before running this script")
    print(" For Ubuntu: ' sudo apt-get install clustalw '")
    quit(save = "no", status = 1, runLast = FALSE)
}

nhmmer.location <- Sys.which("nhmmer")
if(nhmmer.location == ""){
    print("Please install the 'HMMer' software suite before running this script")
    print(" For Ubuntu: ' sudo apt-get install hmmer '")
    quit(save = "no", status = 1, runLast = FALSE)
}




## +++ Inline functions:

align.clustal <- function(mfas, seqtype = NULL, outfile = NULL, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = NULL, nr.cpu = 1) {
	## uwe.menzel@slu.se  uwe.menzel@gmail.com
	clustalw.location <- Sys.which("clustalo")
	if(clustalw.location == "") stop("\n\n  align.clustal: program 'clustalw' not found.\n\n")
	if(!file.exists(mfas)) stop(paste("\n\n  align.clustal: file ' ", mfas, " ' does not exist.\n\n", sep=""))
	if(is.null(outfile)) outfile = paste(mfas, "align.phy", sep="_")
	if(is.null(seqtype)) stop("\n\n  align.clustal: argument 'seqtype' must be one of 'DNA' or 'PROTEIN' \n\n")
	if(!seqtype %in% c("DNA", "PROTEIN")) stop("\n\n  align.clustal: argument 'seqtype' must be one of 'DNA' or 'PROTEIN' \n\n")
	if(!outorder %in% c("ALIGNED", "INPUT")) stop("\n\n  align.clustal: argument 'outorder' must be one of 'ALIGNED' or 'INPUT' \n\n")
	if(!outform %in% c("CLUSTAL", "GCG", "GDE", "PHYLIP", "PIR", "NEXUS", "FASTA"))
		stop("\n\n  align.clustal: argument 'outform' must be one of 'CLUSTAL', 'GCG', 'GDE', 'PHYLIP', 'PIR', 'NEXUS' or 'FASTA' \n\n")
	if(is.null(logf)) logf = paste(gsub(":", "_", gsub(" ", "_", date())), "clustalw.log", sep="_")	# "Wed_Jun_21_11_30_46_2017_clustalw.log"
	#command = paste("clustalw -INFILE=", mfas, " -ALIGN -TREE -TYPE=", seqtype, " -OUTPUT=", outform, " -BOOTSTRAP=", nr.boot, " -OUTFILE=", outfile, "  > ", logf, " 2>&1", sep="")
	command = paste("clustalo --threads=", nr.cpu, " -i ", mfas, " -o ", outform, " > ", logf, " 2>&1", sep="")
	cat(paste("  Running clustalo on", basename(mfas), "; Logfile is", basename(logf), "\n"))
	res <- try(system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL))
	if(res != 0) stop("\n\n  align.clustal: Sorry, no success.\n\n")
	result = list()
	result[["input"]] = mfas
	result[["alignment"]] = outfile
	result[["format"]] = outform
	return(result)
}


build.hmm <- function(alignment = NULL, outhmm = NULL, hmm.name = NULL, summary.file = NULL, new.alignment = NULL, informat = "PHYLIP", nr.cpu = 1) {
	## uwe.menzel@slu.se  uwe.menzel@gmail.com
	hmmbuild.location <- Sys.which("hmmbuild")
	if(hmmbuild.location == "") stop("  build.hmm: program 'hmmbuild' not found.\n\n")
    if(!file.exists(alignment)) stop("  build.hmm: alignment file not found.\n\n")
	if(is.null(outhmm)) outhmm = paste(alignment, "hmm", sep=".")
	if(is.null(hmm.name)) hmm.name = "HMM"
	if(is.null(summary.file))  summary.file = paste(hmm.name, "summary.txt", sep = "_")
	if(is.null(new.alignment)) new.alignment = paste(hmm.name, "HMMer_alignment.phy", sep = "_")
	cat(paste("  Building profile HMM from", basename(alignment), " ... writing to", basename(outhmm), "\n"))
	command = paste("hmmbuild --cpu", nr.cpu, "--dna --fast  -n", hmm.name, "-o", summary.file, "-O", new.alignment, "--informat", informat, outhmm,  alignment)
	res <- try(system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL))
	if(res != 0) stop("\n\n  build.hmm: Sorry, no success.\n\n")
	result = list()
	result[["input"]] = alignment
	result[["output"]] = outhmm
	result[["summary"]] = summary.file
	result[["alignment"]] = new.alignment
	return(result)
}


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


start.html <- function(file, title="Test", overwrite=TRUE) {
  if(file.exists(file)) {
    if(overwrite) {
      unlink(file)
    } else {
      stop(paste("\n  File", file, "already exists, overwritting forbidden by user.\n\n"))
    }
  }
  sink(file)
  cat('<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"')
  cat("\n")
  cat('   "http://www.w3.org/TR/html4/loose.dtd">')
  cat("\n\n")

  cat("<html>", sep="\n")
  cat("<head>", sep="\n")
  cat(paste("<title>", title, "</title>"), sep="\n")
  cat("</head>", sep="\n")
  cat("<body>", sep="\n")
  cat("\n\n")

  cat('<br>')
  cat('<div align="center">')
  cat(paste('<h2><font color="blue">', title, '</font></h2><br>'))
  cat('</div>')
  #cat('<br>')
  #cat('<br>')
  #cat('<br>')
  cat("\n\n")
  sink()
}


set.heading <- function(file, heading, color="blue") {
 if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  cat("\n"); cat("<br>")
  cat(paste('<font color="', color, '" size="+1">', heading, '</font>', sep=""))
  cat("\n")
  sink()
}


line.feed <- function(file, N) {
 if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  breaks = rep('<br>', N)
  cat(breaks)
  sink()
}

insert.text <- function(file, text, color="black") {
  if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  cat("\n")
  cat(paste('<font color="', color,'">', text, '</font>'))
  cat("\n")
  sink()
}


insert.blanks <- function(file, N) {
 if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  blanks = rep('&nbsp;', N)
  cat(blanks)
  sink()
}

stop.html <- function(file) {
  if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  cat("\n\n")
  # cat('<br>')
  cat('mailto: <a href="mailto: uwe.menzel@gmail.com">uwe.menzel@gmail.com</a>')
  cat('<br>');cat('<br>');cat('<br>');cat('<br>');cat('<br>')
  cat("\n\n")
  cat("</body>", sep="\n")
  cat("</html>", sep="\n")
  sink()
}

return.stained <- function(text, color="red")  return(paste('<font color="', color, '">', text, '</font>', sep=""))

insert.line <- function(file, percent=100, size=1) {
 if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  cat("\n"); cat("<br>")
  cat(paste('<hr width="', percent, '%" size="', size, '">', sep=""))
  cat("<br>");cat("\n")
  sink()
}


set.link <- function(file, target, caption="Link") {
 if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  link.tag = paste('<a href="', target, '">', caption, "</a>", sep="")
  cat(link.tag)
  cat("\n")
  sink()
}

put.warning <- function(file, message) {
  if(!file.exists(file)) {
    stop(paste("\n  File", file, "does not exist.\n\n"))
  } else {
    sink(file, append=T)
  }
  cat("\n")
  cat(paste('<font color="red">', paste("  WARNING:", message), '</font><br>'))
  cat("\n")
  sink()
}




## +++ Input parameters:

args <- commandArgs(TRUE)
cat("\n")
if (length(args) < 2) {
  cat("\n  Usage:    Rscript  --vanilla  quantizyme_model.R    <project-ID>   <transcript(fasta)>    [subtree_list (txt)]\n")
  cat("  Example:  Rscript  --vanilla  quantizyme_model.R   GT48    GT48_enzyme.fasta     subtree_list.txt  \n\n")
  cat("    <project-ID>   : some identifier you choose to ease identification of the output.\n")
  cat("    <transcript>   : fasta file containing the DNA sequences you want to match to the metagenomic probe.\n")
  cat("    <subtree_list> : (optional) a text file specifying the subtrees (a separate HMM model is made for each subtree).\n")
  cat("      NOTE: the configuration file ' quant_config.R ' will also be loaded.\n\n\n")
  quit("no")
}

subtree_by_list = ifelse(length(args) == 3, TRUE, FALSE)

project.id = args[1]
transcript.fasta = args[2]
if(subtree_by_list) subtrees.txt = args[3]
cat(paste("  Project ID:", project.id, "\n"))

if(file.exists(transcript.fasta)) {
	cat(paste("  Input fasta file:", transcript.fasta, "\n"))
	transcript.fas = read.fasta(transcript.fasta, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
	if(class(transcript.fas) != "list") stop(paste("\n Reading file '", transcript.fasta, "' not successful. \n\n"))
	nr.seqs = length(transcript.fas)
} else {
	stop(paste(" File '", transcript.fasta, "' not found. \n\n"))
}


if(subtree_by_list) {
	if(file.exists(subtrees.txt)) {
		cat(paste("  Input subtrees file:", subtrees.txt, "\n"))
	} else {
		stop(paste(" File '", subtrees.txt, "' not found. \n\n"))
	}
	subtrees <- read.table(subtrees.txt, stringsAsFactors = FALSE)
	if(ncol(subtrees) != 2) stop(paste("\n Subtree file '", subtrees.txt, "' must have two colums. \n\n"))
	if(class(unique(subtrees[,2])) != "integer") stop(paste("\n Second column of '", subtrees.txt, "' must be integers indicating the subtree. \n\n"))
}

if(subtree_by_list) {
	if(length(transcript.fas) != nrow(subtrees))
		stop(paste("  Input fasta has", length(transcript.fas), "sequences, but '", subtrees.txt, "' has", nrow(subtrees), "entries.\n\n"))
}

if(file.exists("quant_config.R")) {
	cat(paste("  Loading configuration file: ' quant_config.R ' \n"))
	source("quant_config.R")
} else {
	stop(paste(" Configuration file ' quant_config.R ' not found. \n\n"))
}
cat("\n")



# Check if the mandatory variables were defined in the configuration file:
if(exists("nr.boot") && (nr.boot != "")) {
	cat(paste("  Number of bootstrap samples in clustalw:", nr.boot, "\n"))
} else {
	stop(" Mandatory parameter ' nr.boot ' not assigned in the configuration file.\n\n")
}

if(exists("nr.trials.random.picking") && (nr.trials.random.picking != "")) {
	cat(paste("  Number of trials in random picking from subgroups:", nr.trials.random.picking, "\n"))
} else {
	stop(" Mandatory parameter ' nr.trials.random.picking ' not assigned in the configuration file.\n\n")
}

if(exists("subgroup.percent") && (subgroup.percent != "")) {
	cat(paste("  Approximate size of a subgroup in a subtree:", subgroup.percent, "% of the subtree size.\n"))
} else {
	stop(" Mandatory parameter ' subgroup.percent ' not assigned in the configuration file.\n\n")
}

if(exists("hmmbuild.nr.cpu") && (hmmbuild.nr.cpu != "")) {
	cat(paste("  Number of CPU's to be used in hmmbuild:", hmmbuild.nr.cpu, "\n"))
} else {
	stop(" Mandatory parameter ' hmmbuild.nr.cpu ' not assigned in the configuration file.\n\n")
}





## +++ Create output folder:

outfolder = paste(project.id, "MODEL", sep="_")
if(file.exists(outfolder)) {
	if(interactive())  {
		answ <- readline(prompt=paste("  Output folder", outfolder, "already exists. Overwrite?(Y/n): "))
	} else {
		cat(paste("\n  Output folder", outfolder, "already exists. Overwrite?(Y/n): "))
		answ = readLines("stdin", n=1)
	}
	if(tolower(substr(answ,1,1)) == "y") {
		res = unlink(outfolder, recursive = TRUE)
		if(res != 0) stop(paste("  Could not remove folder '", outfolder, "'\n\n"))
		res = dir.create(file.path(getwd(), outfolder), showWarnings = FALSE)		# create the directory
		if(!res) stop(paste("  Could not create folder '", outfolder, "'\n\n"))
	} else {
		stop(" Script stopped on user request to prevent overwriting of existing folder.\n\n")
	}
} else {
	dir.create(file.path(getwd(), outfolder), showWarnings = FALSE)		# create the directory
}






## +++ Initialize output HTML:

outhtml = paste(project.id, "model.html", sep="_")
outhtml = file.path(getwd(), outfolder, outhtml)
startfolder = getwd()
input.fasta = file.path(getwd(), transcript.fasta)
cat(paste("  Results will be linked in '", outhtml, "'\n"))
cat("     (open this file in a web browser).\n\n")

start.html(outhtml, title="Quantify DNA in metagenomic probes", overwrite=T)
set.heading(outhtml, heading = "Model definition", color="blue")
set.heading(outhtml, heading = paste("Project ID:", project.id), color="blue")
line.feed(outhtml, 1)
line.feed(outhtml, 1)

insert.text(outhtml, text=paste(" Start folder:", getwd()), color="black")
line.feed(outhtml, 1)
insert.text(outhtml, text=date(), color="black")
line.feed(outhtml, 1)
insert.text(outhtml, text=R.Version()$version.string, color="black")
line.feed(outhtml, 1)
insert.line(outhtml, percent=100, size=1)

prog.name = "quantizyme_model.R"
insert.text(outhtml, text=paste(" Script executed:", return.stained(prog.name, color="red")), color="black")
line.feed(outhtml, 1)
insert.text(outhtml, text=paste(" Command line arguments:", paste(args, collapse="&nbsp;&nbsp; &nbsp;"), sep = "&nbsp;&nbsp; &nbsp;"), color="black")
line.feed(outhtml, 1)

insert.text(outhtml, text="Configuration file:", color="black")
res <- file.copy("quant_config.R", file.path(getwd(), outfolder))
config.file = file.path(getwd(), outfolder, "quant_config.R")
set.link(outhtml, config.file, caption="quant_config.R")
line.feed(outhtml, 1)

insert.text(outhtml, text="Input fasta: ", color="black")
set.link(outhtml, input.fasta, caption=transcript.fasta)
insert.blanks(outhtml, 1)
insert.text(outhtml, text=paste("(", return.stained(length(transcript.fas), "blue"), "sequences )"), color="black")
line.feed(outhtml, 1)

insert.text(outhtml, text=paste(" Output folder: ", return.stained(outfolder, color="red")), color="black")
line.feed(outhtml, 1)
insert.line(outhtml, percent=100, size=1)


# Parameters read from quant_config.R
insert.text(outhtml, text="Parameters read from ' quant_config.R '", color="black")
line.feed(outhtml, 1)

insert.blanks(outhtml, 2)
insert.text(outhtml, text = paste("  Number of bootstrap samples in clustalw:", return.stained(nr.boot, color="blue")), color="black")
line.feed(outhtml, 1)

insert.blanks(outhtml, 2)
insert.text(outhtml, text = paste("  Number of trials in random picking from subgroups:", return.stained(nr.trials.random.picking, color="blue")), color="black")
line.feed(outhtml, 1)

insert.blanks(outhtml, 2)
insert.text(outhtml, text = paste("  Approximate size of a subgroup in a subtree:", return.stained(subgroup.percent, color="blue"), "% of the subtree size."), color="black")
line.feed(outhtml, 1)

insert.blanks(outhtml, 2)
insert.text(outhtml, text = paste("  Number of CPU's used when creating the profile HMM (hmmbuild command):", return.stained(hmmbuild.nr.cpu, color="blue")), color="black")
line.feed(outhtml, 1)
insert.line(outhtml, percent=100, size=1)

insert.blanks(outhtml, 2)
insert.text(outhtml, text = paste("  Number of CPU's used when creating the MSA (clustalO command):", return.stained(clustalo.nr.cpu, color="blue")), color="black")
line.feed(outhtml, 1)
insert.line(outhtml, percent=100, size=1)


## +++ Remove very short and very long reads from the input file:

seqlength = as.integer(lapply(transcript.fas, nchar))
names(seqlength) = names(transcript.fas)
nr.seqs = length(transcript.fas)

# Display sequence length distribution on screen
X11(width=6, height=6)
hist(seqlength, col="red", breaks=25, main="Transcript sequence length distribution", font.main=1, xlab="Length")
mtext(side=3, paste("Based on", nr.seqs, "sequences"), col="blue", cex=0.9)

# save plot to file
fn = paste(project.id, "read_length_distribution.pdf", sep="_")
#cat(paste("  Saving histogram for transcript sequence length to", fn, "\n"))
outpdf = file.path(getwd(), outfolder, fn)
#cat(paste(" Opening device", "\n"))
#png(file=outpng)
#cat(paste(" Plotting", "\n"))
#hist(seqlength, col="red", breaks=25, main="Transcript sequence length distribution", font.main=1, xlab="Length")
#mtext(side=3, paste("Based on", nr.seqs, "sequences"), col="blue", cex=0.9)
#dev.off()

### try with savePlot
#savePlot(outpng, type = "png")
#####

### original
invisible(dev.copy2pdf(file = outpdf, out.type = "pdf"))
#invisible(dev.off())
#####

# cat(paste("  Opened devices (1)", fn, "\n"))
# dev.list()
# cat(paste("  dev copy", fn, "\n"))
##dev.copy(png, file = outpng, unit="px", width=960, height=960, res=120)
# cat(paste("  Opened devices (2)", fn, "\n"))
# dev.list()
##dev.off()

cat(paste("  Histogram for transcript sequence length saved to", fn, "\n"))

insert.text(outhtml, text="Fasta read length distribution:", color="black")
set.link(outhtml, outpdf, caption="histogram")
#set.link(outhtml, outpng, caption="histogram")
line.feed(outhtml, 1)


if(interactive())  {
	answ = readline(prompt = "  Remove short/long sequences?(Y/n): ")
} else {
	cat("  Remove short/long sequences?(Y/n): ")
	answ = readLines("stdin", n=1)
}

if(tolower(substr(answ,1,1)) == "y") {
	repeat {
		if(interactive())  {
			low.thres = as.integer(readline(prompt="  Lower threshold(integer): "))
		} else {
			cat("  Lower threshold(integer): ")
			low.thres = readLines("stdin", n=1)
			low.thres = as.integer(low.thres)
		}

		if(interactive())  {
			up.thres = as.integer(readline(prompt="  Upper threshold(integer): "))
		} else {
			cat("  Upper threshold(integer): ")
			up.thres = readLines("stdin", n=1)
			up.thres = as.integer(up.thres)
		}

		X11(width=6, height=6)
		hist(seqlength, col="red", breaks=25, main="Transcript sequence length distribution", font.main=1, xlab="Length/nt")
		mtext(side=3, " Including chosen thresholds", col="blue", cex=0.9)
		abline(v=low.thres, col="blue", lty=3)
		abline(v=up.thres, col="blue", lty=3)

		nr.skipped.low = sum(seqlength < low.thres)
		nr.skipped.up = sum(seqlength > up.thres)
		cat(paste("  This removes", nr.skipped.low, "sequences at the lower tail and", nr.skipped.up, "sequences at the upper tail.\n"))

		if(interactive())  {
			answ <- readline(prompt="  Threshold accepted?(Y/n): ")
		} else {
			cat("  Threshold accepted?(Y/n): ")
			answ = readLines("stdin", n=1)
		}

		if(tolower(substr(answ,1,1)) == "y") break
	}

	X11(width=6, height=6)
	hist(seqlength, col="red", breaks=25, main="Transcript sequence length distribution", font.main=1, xlab="Length/nt")
	mtext(side=3, paste("Based on", nr.seqs, "sequences"), col="blue", cex=0.9)
	abline(v=low.thres, col="blue", lty=3)
	abline(v=up.thres, col="blue", lty=3)
	# fn = paste(project.id, "thresholded_distribution.png", sep="_")
	# outpng1 = file.path(getwd(), outfolder, fn)
	# invisible(dev.copy(png, file = outpng1, unit="px", width=960, height=960, res=120)); invisible(dev.off())
	fn = paste(project.id, "thresholded_distribution.pdf", sep="_")
	outpdf1 = file.path(getwd(), outfolder, fn)
	invisible(dev.copy2pdf(file = outpdf1)); invisible(dev.off())
	cat(paste("  Histogram for thresholded sequence length saved to", fn, "\n"))

	insert.text(outhtml, text="Length distribution with chosen thresholds:", color="black")
	set.link(outhtml, outpdf1, caption="histogram")
	line.feed(outhtml, 1)
	insert.text(outhtml, text=paste("Lower threshold:", return.stained(low.thres, "blue")), color="black")
	insert.blanks(outhtml, 2)
	insert.text(outhtml, text=paste("Upper threshold:", return.stained(up.thres, "blue")), color="black")
	line.feed(outhtml, 1)

	# extract sequences clipped at lower tail:
	ind.low = which(seqlength < low.thres)
	if(length(ind.low) == 0) {		# nothing to clip
		cat(paste("  No sequences clipped at lower tail.\n"))
	} else {
		to_extract = transcript.fas[ind.low]
		if(length(to_extract) != nr.skipped.low) stop(" Wrong number of skipped sequences at lower tail.\n\n")
		fn = paste(project.id, "skipped_low.fasta", sep = "_")
		outtxt1 = file.path(getwd(), outfolder, fn)
		write.fasta(to_extract, names = names(to_extract), file.out = outtxt1, open = "w", nbchar = 60, as.string = TRUE)
		cat(paste("  Sequences clipped on lower tail written to '", fn, "' \n"))
	}

	# extract sequences clipped at upper tail:
	ind.up = which(seqlength > up.thres)
	if(length(ind.up) == 0) {
		cat(paste("  No sequences clipped at upper tail.\n"))
	} else {
		to_extract = transcript.fas[ind.up]
		if(length(to_extract) != nr.skipped.up) stop(" Wrong number of skipped sequences at upper tail.\n\n")
		fn = paste(project.id, "skipped_up.fasta", sep = "_")
		outtxt2 = file.path(getwd(), outfolder, fn)
		write.fasta(to_extract, names = names(to_extract), file.out = outtxt2, open = "w", nbchar = 60, as.string = TRUE)
		cat(paste("  Sequences clipped at upper tail written to '", fn, "' \n"))
	}

	# Link clipped parts:
	if(length(ind.low) == 0) {
		insert.text(outhtml, text="No sequence clipped at lower tail.", color="black")
		line.feed(outhtml, 1)
	} else {
		insert.text(outhtml, text="Sequence clipped at lower tail:", color="black")
		set.link(outhtml, outtxt1, caption="fasta")
		insert.blanks(outhtml, 1)
		insert.text(outhtml, text=paste("(", return.stained(nr.skipped.low, "blue"), "sequences )"), color="black")
		line.feed(outhtml, 1)
	}

	if(length(ind.up) == 0) {
		insert.text(outhtml, text="No sequence clipped at upper tail.", color="black")
		line.feed(outhtml, 1)
	} else {
		insert.text(outhtml, text="Sequence clipped at upper tail:", color="black")
		set.link(outhtml, outtxt2, caption="fasta")
		insert.blanks(outhtml, 1)
		insert.text(outhtml, text=paste("(", return.stained(nr.skipped.up, "blue"), "sequences )"), color="black")
		line.feed(outhtml, 1)
	}

	ind.to.remove = integer(0)
	if(length(ind.up)  != 0) ind.to.remove = ind.up
	if(length(ind.low) != 0) ind.to.remove = c(ind.to.remove, ind.low)
	if(length(ind.to.remove) != 0) transcript.fas = transcript.fas[-ind.to.remove]

	start.fasta = paste(project.id, "start.fasta", sep = "_")
	start.fasta.location = file.path(getwd(), outfolder, start.fasta)
	write.fasta(transcript.fas, names = names(transcript.fas), file.out = start.fasta.location, open = "w", nbchar = 60, as.string = TRUE)
	cat(paste("  Thresholded fasta written to '", start.fasta, "' \n"))

	insert.text(outhtml, text="Transcript sequence after thresholding:", color="black")
	set.link(outhtml, start.fasta.location, caption="fasta")
	insert.blanks(outhtml, 1)
	insert.text(outhtml, text=paste("(", return.stained(length(transcript.fas), "blue"), "sequences )"), color="black")
	line.feed(outhtml, 1)
} else {
	start.fasta = paste(project.id, "start.fasta", sep = "_")
	start.fasta.location = file.path(getwd(), outfolder, start.fasta)
	success <- file.copy(input.fasta, start.fasta.location, overwrite = TRUE)  # just copy the original fasta to the start.fasta.location
	if(!success) stop(paste(" File '", start.fasta, "' could not be copied.\n\n"))
	insert.blanks(outhtml, 1)
	insert.text(outhtml, text="No length-specific clipping was done (users choice).", color="black")
	line.feed(outhtml, 1)
}

nr.seqs = length(transcript.fas)
cat(paste("\n  Starting analysis with", nr.seqs, "sequences.\n"))
if(!file.exists(start.fasta.location)) stop(" Start file not available.\n\n")



## +++ SUBTREEing
insert.line(outhtml, percent=100, size=1)

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
	insert.text(outhtml, text="Subtreeing is conducted by clustering.", color="black")
	line.feed(outhtml, 1)

	if(exists("phylip.all") && (phylip.all != "")) {
		if(file.exists(phylip.all)) {
			cat(paste("  Alignment for all sequences will be loaded from '", phylip.all, "' (specified in configuration file)\n"))
			aln <- read.alignment(phylip.all, format="phylip", forceToLower = FALSE)
			if(class(aln) != "alignment") stop(paste(" Could not read alignment '", out.phy, "'.\n\n"))
			if(aln$nb != length(transcript.fas))
				stop(paste(" Alignment has", aln$nb, "entries, but transcript had", length(transcript.fas), "entries.\n\t Did you specify the correct alignment file?\n\n"))
			success <- file.copy(phylip.all, file.path(getwd(), outfolder), overwrite = TRUE)
			if(!success) stop(paste(" File '", phylip.all, "' could not be copied to folder", outfolder, ".\n\n"))
			insert.text(outhtml, text="Loaded alignment for transcript sequence:", color="black")
			set.link(outhtml, file.path(getwd(), outfolder, phylip.all), caption="phylip format")
			line.feed(outhtml, 1)
		} else {
			stop(paste(" Alignment file '", phylip.all, "' was specified in configuration file, but was not found in the current folder.\n\n"))
		}
	} else {
		fn = paste(project.id, "align1.phy", sep="_")
		out.phy = file.path(getwd(), outfolder, fn)
		clustal.log = paste(project.id, "align1_clustal.log", sep="_")
		fn = file.path(getwd(), outfolder, clustal.log)
		aln.res <- align.clustal(start.fasta.location, seqtype = "DNA", outfile = out.phy, outform = "PHYLIP", outorder = "ALIGNED", nr.boot = 100, logf = fn)
		if(!file.exists(out.phy))  stop(" Could not complete clustalw alignment.\n\n")
		aln <- read.alignment(out.phy, format="phylip", forceToLower = FALSE)
		if(class(aln) != "alignment") stop(paste(" Could not read alignment '", out.phy, "'.\n\n"))
		if(aln$nb != length(transcript.fas)) stop(paste(" Alignment has", aln$nb, "entries, but transcript had", length(transcript.fas), "entries.\n\n"))
		insert.text(outhtml, text="Calculated alignment for transcript sequence:", color="black")
		set.link(outhtml, out.phy, caption="phylip format")
		line.feed(outhtml, 1)
	}

	aln.1 = as.phyDat(aln)
	if(class(aln.1) != "phyDat") stop(" Could not convert alignment to phyDat format.\n\n")

	distmat = dist.ml(aln.1, model="JC69")
	if(class(distmat) != "dist") stop(" Could not calculate distance matrix from alignment.\n\n")
	if(sum(is.na(distmat))	!= 0) stop(" Distance matrix includes missing data.\n\n")

	cluster = hclust(distmat, method = "ward.D")
	if(class(cluster) != "hclust") stop(" Hierarchical clustering failed.\n\n")

	if(length(cluster$labels) < 25) {
		labels = gsub(" ", "", cluster$labels)
	} else {
		labels = FALSE
	}
	if(is.logical(labels)) cat("  Labels in the phylogenetic trees cannot be shown: too many sequences.\n\n")

	X11(width=6, height=6)
	plot(cluster, labels = labels, hang = 0.1, ann = TRUE, main = paste("Phylogenetic tree, project", project.id), sub = "", xlab = "", ylab = "", font.main=1)
	mtext(side=1, paste(length(transcript.fas), "sequences"), col="blue", cex=0.9)
	# fn = paste(project.id, "hc_tree_1.png", sep="_")
	# outpng3 = file.path(getwd(), outfolder, fn)
	# invisible(dev.copy(png, file = outpng3, unit="px", width=960, height=960, res=120)); invisible(dev.off())
	fn = paste(project.id, "hc_tree_1.pdf", sep="_")
	outpdf3 = file.path(getwd(), outfolder, fn)
	invisible(dev.copy2pdf(file = outpdf3, out.type = "pdf")); invisible(dev.off())
	cat(paste("  Phylogenetic tree (version 1) saved to", fn, "\n"))
	# link plot
	insert.text(outhtml, text="Phylogenetic tree, rectangle format (hclust):", color="black")
	set.link(outhtml, outpdf3, caption="pdf")
	insert.blanks(outhtml, 2)
	# link tree in Newick format:
	cluster.phylo <- as.phylo(cluster)
	fn = "tree1.new"
	tree.file = file.path(getwd(), outfolder, fn)
    write.tree(phy = cluster.phylo, file = tree.file)
    cat(paste("  Neighbor joining tree (version 1) saved to '", fn, "'\n\n"))
	set.link(outhtml, tree.file, caption="Newick")
	line.feed(outhtml, 1)

    cluster2 = nj(distmat)
	if(class(cluster2) != "phylo") stop(" Could not calculate phylogenetic tree by neighbor-joining.\n\n")

    if(cluster2$Nnode < 25) {
		labels = TRUE
	} else {
		labels = FALSE
	}

	X11(width=6, height=6)
	plot.phylo(cluster2, type = "unrooted", show.tip.label = labels, main = paste("Phylogenetic tree, project", project.id), font.main = 1)
	fn = paste(project.id, "hc_tree_2.png", sep="_")
	outpdf7 = file.path(getwd(), outfolder, fn)
	invisible(dev.copy2pdf(file = outpdf7, out.type = "pdf")); invisible(dev.off())
	cat(paste("  Phylogenetic tree (version 2) saved to", fn, "\n"))
	# link plot
	insert.text(outhtml, text="Phylogenetic tree, radial format (nj, ape):", color="black")
	set.link(outhtml, outpdf7, caption="pdf")
	insert.blanks(outhtml, 2)
	# link tree in Newick format:
	fn = "tree2.new"
	tree.file = file.path(getwd(), outfolder, fn)
    write.tree(phy = cluster2, file = tree.file)
    cat(paste("  Neighbor joining tree (version 2) saved to '", fn, "'\n\n"))
    set.link(outhtml, tree.file, caption="Newick")
	line.feed(outhtml, 1)

	options(warn=-1)
	repeat {
		if(interactive())  {
			nrHMM <- readline(prompt="  Number of Hidden Markov Models to create?(integer): ")
		} else {
			cat("  Number of Hidden Markov Models to create?(integer): ")
			nrHMM = readLines("stdin", n=1)
		}
		nrHMM = as.integer(nrHMM)
		if(!is.na(nrHMM))  {
			if(nrHMM == round(nrHMM)) break
		}
	}
	options(warn=0)

	cat(paste("\n  Building", nrHMM, "profile Hidden Markov Models:\n"))
	insert.line(outhtml, percent=100, size=1)
	insert.text(outhtml, text=paste("Number of HMM models to build:", return.stained(nrHMM, "red")), color="black")
	line.feed(outhtml, 1)

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

	# link the results on website:
	for (i in 1:nrHMM) {
		insert.blanks(outhtml, 2)
		insert.text(outhtml, text=paste("  Subtree", i, ":"), color="black")
		fn = file.path(getwd(), outfolder, subtree.fn[[i]])
		set.link(outhtml, fn, caption="fasta")
		insert.blanks(outhtml, 2)
		insert.text(outhtml, text=paste("(", length(subtree[[i]]), "sequences" , ")"), color="black")
		line.feed(outhtml, 1)
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
	insert.line(outhtml, percent=100, size=1)
	warn.txt = paste("  The following subtrees have been removed because they contain less than 2 sequences:", paste(" Subtree", ind.to.remove, collapse=" , "))
	put.warning(outhtml, warn.txt)
	if(nrHMM == 0) stop("  No subtrees left. Exiting.\n\n")
}



## +++ Create multiple alignment for each of the subtrees:

insert.line(outhtml, percent=100, size=1)
insert.text(outhtml, text=paste("Calculation of the alignments for the", nrHMM, "subtrees:"), color="black")
line.feed(outhtml, 1)

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
	insert.blanks(outhtml, 2)
	insert.text(outhtml, text=paste("Calculated alignment for subtree", i, ":"), color="black")
	set.link(outhtml, out.phy, caption="phylip format")
	line.feed(outhtml, 1)
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

insert.line(outhtml, percent=100, size=1)
insert.text(outhtml, text="Dendograms for the subtrees:", color="black")
line.feed(outhtml, 1)

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
	insert.blanks(outhtml, 2)
	insert.text(outhtml, text=paste("Phylogenetic tree for subtree", i, ": "), color="black")
	set.link(outhtml, fn, caption="png")
	line.feed(outhtml, 1)
}


## +++ Cut subtree into appropiate number of subgroups:

insert.line(outhtml, percent=100, size=1)
insert.text(outhtml, text="Picking from subtrees ('Thinning'):", color="black")
line.feed(outhtml, 1)

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
		insert.blanks(outhtml,4)
		insert.text(outhtml, text=paste("Picked from subtree", i, ":"), color="black")
		set.link(outhtml, fn, caption="fas")
		insert.blanks(outhtml,2)
		insert.text(outhtml, text=paste("( no trials made - too few sequences )"), color="black")
		line.feed(outhtml, 1)
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
		insert.blanks(outhtml, 2)
		insert.text(outhtml, text=paste("Number of sequences in the groups of subtree", i, ": "), color="black")
		set.link(outhtml, fn, caption="png")
		line.feed(outhtml, 1)
		to_pick = matrix(NA, nrow = num.clust, ncol = nr.trials.random.picking)
		insert.blanks(outhtml,4)
		insert.text(outhtml, text=paste("Picked from subtree", i, ":"), color="black")
		for (k in 1:nr.trials.random.picking) {
			to_pick[,k] = pick_one_from_each_cluster(subtree.group[[i]], silent = TRUE)
			fasta.picked = subtree[[i]][to_pick[,k]]
			if(length(fasta.picked) != num.clust) stop(paste(" Fasta incorrectly picked: i =",i, ", k =", k,"\n\n"))
			picked.fn[i,k] = paste("subtree", i, "trial", k, "picked.fasta", sep="_")
			fn = file.path(getwd(), outfolder, picked.fn[i,k])
			write.fasta(fasta.picked, names = names(fasta.picked), file.out = fn, open = "w", nbchar = 60, as.string = TRUE)
			# link on webpage:
			insert.blanks(outhtml,1)
			insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
			set.link(outhtml, fn, caption="fas")
			if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
		}
		colnames(to_pick) = paste("trial", 1:nr.trials.random.picking, sep="_")
		rownames(to_pick) = paste("cl", 1:num.clust, sep="_")
		if(sum(is.na(to_pick)) != 0) stop("  Could not find indices for sequences to pick from subgroups.\n\n")
	}
}
cat("\n")


## Create profile HMM from each of the fasta files:

insert.line(outhtml, percent=100, size=1)
insert.text(outhtml, text="Create profile Hidden Markov Models for subtrees:", color="black")
line.feed(outhtml, 1)

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
		insert.blanks(outhtml,2)
		insert.text(outhtml, text=paste("HMM for subtree", i, ":"), color="black")
		set.link(outhtml, out.hmm, caption="txt")
		insert.blanks(outhtml,2)
		insert.text(outhtml, text=paste("( no trials made - too few sequences )"), color="black")
		line.feed(outhtml, 1)
	} else {
		insert.blanks(outhtml,2)
		insert.text(outhtml, text=paste("HMM's for subtree", i, ":"), color="black")
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
			insert.blanks(outhtml,1)
			insert.text(outhtml, text=paste("trial ", k, ":", sep=""), color="black")
			set.link(outhtml, out.hmm, caption="txt")
			if(k == nr.trials.random.picking) line.feed(outhtml, 1) else insert.blanks(outhtml, 1)
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
insert.line(outhtml, percent=100, size=1)
insert.text(outhtml, text = "Hidden Markov Models:", color = "black")
zip.location = file.path(getwd(), outfolder, zipfile)
line.feed(outhtml, 1)
insert.blanks(outhtml,2)
insert.text(outhtml, text="Models:", color="black")
insert.blanks(outhtml,1)
set.link(outhtml, zip.location, caption = zipfile)
insert.text(outhtml, text=paste("(", return.stained(number.models, "blue"), "models )"), color="black")
line.feed(outhtml, 1)
insert.blanks(outhtml,2)
insert.text(outhtml, text="Listing of models + parameter file:", color="black")
insert.blanks(outhtml,1)
list.location = file.path(getwd(), outfolder, zip.list)
set.link(outhtml, list.location, caption = zip.list)
insert.blanks(outhtml, 1)
line.feed(outhtml, 1)


## +++ End
insert.line(outhtml, percent=100, size=1)
stop.html(outhtml)
cat(paste("\n  Results written to '", outhtml, "'\n\n"))

# uwe.menzel@matstat.de
# uwe.menzel@gmail.com
