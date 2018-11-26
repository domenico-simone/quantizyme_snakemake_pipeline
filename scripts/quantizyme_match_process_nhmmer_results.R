required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs) 
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))  		
suppressMessages(library(phangorn))		
suppressMessages(library(ape))			 

search.hmm <- function(hmm = NULL, fasta = NULL, out.tbl = NULL, nr.cpu = 1, eval = 1, alignments = NULL, logf = NULL) {
	## uwe.menzel@slu.se  uwe.menzel@gmail.com  
	nhmmer.location <- Sys.which("nhmmer")
	if(nhmmer.location == "") stop("  search.hmm: program 'nhmmer' not found.\n\n")
	if(is.null(hmm)) stop("  search.hmm: ' hmm ' is a mandatory parameter.\n\n")
	if(is.null(fasta)) stop("  search.hmm: ' fasta ' is a mandatory parameter.\n\n")
	if(is.null(out.tbl)) stop("  search.hmm: ' out.tbl ' is a mandatory parameter.\n\n")
	if(is.null(alignments))  alignments = "nhmmer_alignment.msa" 
	if(is.null(logf)) logf = paste(gsub(":", "_", gsub(" ", "_", date())), "nhmmer.log", sep="_")
	cat(paste("  Searching", basename(fasta), "with", basename(hmm), "\n"))
    command = paste("nhmmer --cpu", nr.cpu, "-o", logf, "-A", alignments, "--tblout", out.tbl, "-E", eval, "--tformat fasta", hmm, fasta)	
	res <- try(system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL))
	if(res != 0) stop("\n\n  search.hmm: Sorry, no success.\n\n")
	result = list()
	result[["hmm"]] = hmm
	result[["fasta"]] = fasta
	result[["table"]] = out.tbl
	result[["eval"]] = eval
	result[["alignments"]] = alignments
	return(result)
} 