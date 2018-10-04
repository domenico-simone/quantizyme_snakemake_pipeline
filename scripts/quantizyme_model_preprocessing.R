
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
