required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs)
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))
suppressMessages(library(phangorn))
suppressMessages(library(ape))

args <- commandArgs(TRUE)
#project.id = args[1]
outplot = args[1]
transcript.fasta = args[2]
#if(subtree_by_list) subtrees.txt = args[3]
#cat(paste("  Project ID:", project.id, "\n"))

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

# save plot to file
# outfolder is {reference_transcripts_length_distribution}/{projectID}_transcript_length_distribution.pdf
###outfolder = "reference_transcripts_length_distribution"
#fn = paste(project.id, "transcript_length_distribution.pdf", sep="_")
###fn = file.path(outfolder, paste(project.id, "transcript_length_distribution.pdf", sep="_"))
#cat(paste("  Saving histogram for transcript sequence length to", fn, "\n"))
#outpdf = file.path(getwd(), outfolder, fn)
outpdf = outplot
pdf(outpdf)
hist(seqlength, col="red", breaks=25, main="Transcript sequence length distribution", font.main=1, xlab="Length")
mtext(side=3, paste("Based on", nr.seqs, "sequences"), col="blue", cex=0.9)
dev.off()
