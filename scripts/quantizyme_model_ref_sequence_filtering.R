required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs)
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))
suppressMessages(library(phangorn))
suppressMessages(library(ape))
suppressMessages(library(optparse))

# args <- commandArgs(TRUE)
# project.id = args[1]
# transcript.fasta = args[2]
# low.thres = as.integer(args[3])
# up.thres = as.integer(args[4])

# mandatory opts:
# - projectID
# - input
option_list = list(
    make_option(c("-p", "--projectID"), type="character", default=NULL,
            help="projectID [default= %default] (MANDATORY)", metavar="character"),
	make_option(c("-i", "--input"), type="character", default=NULL,
            help="input file name [default= %default] (MANDATORY)", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL,
            help="output file name [default= %default] (MANDATORY)", metavar="character"),
    make_option(c("-f", "--filter_seqs"), default=FALSE,
            help="filter short/long transcript [default= %default]"),
    make_option(c("-l", "--filter_lower_t"), type="integer", default=0,
            help="lower threshold for transcript filtering if filter_seqs == TRUE"),
    make_option(c("-u", "--filter_upper_t"), type="integer", default=0,
            help="upper threshold for transcript filtering if filter_seqs == TRUE"),
    make_option("--outplot", type="character", default=NULL,
            help="plot file of sequence length thresholded distribution if filter_seqs == TRUE"),
    make_option(c("-d", "--outdir"), type="character", default=NULL,
            help="outdir"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#### Check mandatory parameters
if (is.null(opt$projectID)){
  print_help(opt_parser)
  stop("Please provide projectID.\n", call.=FALSE)
}

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Please provide input file name.\n", call.=FALSE)
}

if (is.null(opt$output)){
  print_help(opt_parser)
  stop("Please provide output file name.\n", call.=FALSE)
}

if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Please provide outdir name.\n", call.=FALSE)
}

if (isTRUE(opt$filter_seqs)){
    if (opt$filter_lower_t == 0 & opt$filter_upper_t == 0){
        print_help(opt_parser)
        stop("Length cutoffs for transcript filtering are not valid.\n", .call=FALSE)
    }
}

#### Assign variables

transcript.fasta = opt$input
project.id = opt$projectID
start.fasta = opt$output # output with filtered sequences
outplot = opt$outplot
low.thres = opt$filter_lower_t
up.thres = opt$filter_upper_t
outdir = opt$outdir

# If no filtering is set, just copy the transcript file to the output file
if (!isTRUE(opt$filter_seqs)){
    file.copy(transcript.fasta, start.fasta)
    stop("", .call=TRUE)
}

# Read input fasta file
if(file.exists(transcript.fasta)) {
	cat(paste("  Input fasta file:", transcript.fasta, "\n"))
	transcript.fas = read.fasta(transcript.fasta, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)
	if(class(transcript.fas) != "list") stop(paste("\n Reading file '", transcript.fasta, "' not successful. \n\n"))
	nr.seqs = length(transcript.fas)
} else {
	stop(paste(" File '", transcript.fasta, "' not found. \n\n"))
}

# get sequence lengths
seqlength = as.integer(lapply(transcript.fas, nchar))
names(seqlength) = names(transcript.fas)
nr.seqs = length(transcript.fas)

# make a new histogram of length distribution only if filtering is enabled
if (isTRUE(opt$filter_seqs)){
    outpdf = outplot
    pdf(outpdf)
    hist(seqlength, col="red", breaks=25, main="Transcript sequence length distribution", font.main=1, xlab="Length/nt")
    mtext(side=3, paste("Based on", nr.seqs, "sequences"), col="blue", cex=0.9)
    abline(v=low.thres, col="blue", lty=3)
    abline(v=up.thres, col="blue", lty=3)
    dev.off()
}

######################################################

nr.skipped.low = sum(seqlength < low.thres)
nr.skipped.up = sum(seqlength > up.thres)
cat(paste("  This removes", nr.skipped.low, "sequences at the lower tail and", nr.skipped.up, "sequences at the upper tail.\n"))

# extract sequences clipped at lower tail:
ind.low = which(seqlength < low.thres)
if(length(ind.low) == 0) {		# nothing to clip
    cat(paste("  No sequences clipped at lower tail.\n"))
} else {
    to_extract = transcript.fas[ind.low]
    if(length(to_extract) != nr.skipped.low) stop(" Wrong number of skipped sequences at lower tail.\n\n")
    fn = paste(project.id, "skipped_low.fasta", sep = "_")
    #outtxt1 = file.path(getwd(), outfolder, fn)
    outtxt1 = file.path(outdir, fn)
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
    outtxt2 = file.path(outdir, fn)
    #outtxt2 = file.path(getwd(), outfolder, fn)
    write.fasta(to_extract, names = names(to_extract), file.out = outtxt2, open = "w", nbchar = 60, as.string = TRUE)
    cat(paste("  Sequences clipped at upper tail written to '", fn, "' \n"))
}

# Link clipped parts:
if(length(ind.low) == 0) {
    cat("No sequence clipped at lower tail.")
    # insert.text(outhtml, text="No sequence clipped at lower tail.", color="black")
    # line.feed(outhtml, 1)
} else {
    # insert.text(outhtml, text="Sequence clipped at lower tail:", color="black")
    cat("Sequences clipped at lower tail.")
    # set.link(outhtml, outtxt1, caption="fasta")
    # insert.blanks(outhtml, 1)
    # insert.text(outhtml, text=paste("(", return.stained(nr.skipped.low, "blue"), "sequences )"), color="black")
    # line.feed(outhtml, 1)
}

if(length(ind.up) == 0) {
    cat("No sequence clipped at upper tail.")
    # insert.text(outhtml, text="No sequence clipped at upper tail.", color="black")
    # line.feed(outhtml, 1)
} else {
    cat("Sequences clipped at upper tail.")
    # insert.text(outhtml, text="Sequence clipped at upper tail:", color="black")
    # set.link(outhtml, outtxt2, caption="fasta")
    # insert.blanks(outhtml, 1)
    # insert.text(outhtml, text=paste("(", return.stained(nr.skipped.up, "blue"), "sequences )"), color="black")
    # line.feed(outhtml, 1)
}

ind.to.remove = integer(0)
if(length(ind.up)  != 0) ind.to.remove = ind.up
if(length(ind.low) != 0) ind.to.remove = c(ind.to.remove, ind.low)
if(length(ind.to.remove) != 0) transcript.fas = transcript.fas[-ind.to.remove]

#start.fasta = paste(project.id, "start.fasta", sep = "_")
#start.fasta.location = file.path(getwd(), outfolder, start.fasta)
write.fasta(transcript.fas, names = names(transcript.fas), file.out = start.fasta, open = "w", nbchar = 60, as.string = TRUE)
cat(paste("  Thresholded fasta written to '", start.fasta, "' \n"))

# insert.text(outhtml, text="Transcript sequence after thresholding:", color="black")
# set.link(outhtml, start.fasta.location, caption="fasta")
# insert.blanks(outhtml, 1)
# insert.text(outhtml, text=paste("(", return.stained(length(transcript.fas), "blue"), "sequences )"), color="black")
# line.feed(outhtml, 1)
