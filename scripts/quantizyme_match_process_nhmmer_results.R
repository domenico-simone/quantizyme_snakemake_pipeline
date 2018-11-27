required.libs = c("seqinr", "phangorn", "ape")
installed.libs = rownames(installed.packages())
missing.libs = setdiff(required.libs, installed.libs) 
if(length(missing.libs) > 0) stop(paste("  Missing the following R-libraries:", paste(missing.libs, collapse = "  "), "\n\n"))

suppressMessages(library(seqinr))          
suppressMessages(library(phangorn))        
suppressMessages(library(ape))             
suppressMessages(library(optparse))

# Quite surely we don't need this function anymore
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


parse_nhmmer_matchtable <- function(match.table, match.fn = NULL, eval = 1) {
    ## uwe.menzel@gmail.com
    if(!file.exists(match.table)) stop(paste("  parse_nhmmer_matchtable: File '", match.table, "' not found.\n\n"))  
    trial = try(read.table(match.table, stringsAsFactors = FALSE), silent = TRUE)    
    if(inherits(trial, "try-error")) {
        cat(paste("  Number of unique matching reads with eval <=", eval, "in table '", basename(match.table), "' : ", 0, "\n"))
        result = list()
        result[["matchtable"]] = match.table
        result[["nr.matches"]] = 0
        result[["matched_reads"]] = NA
        return(result)
    } else {
        matches = trial    # 
        if(class(matches) != "data.frame") stop(paste("  parse_nhmmer_matchtable: File '", match.table, "' not correctly read.\n\n"))
    }
    
    ev = matches[,13]        
    if(!is.numeric(ev)) stop(paste("  parse_nhmmer_matchtable: File '", match.table, "' not correctly parsed (e-values).\n\n"))
    matches = matches[which(ev <= eval),]      
    if(nrow(matches) == 0) {
        cat(paste("  Number of unique matching reads with eval <=", eval, "in table '", basename(match.table), "' : ", 0, "\n"))
        result = list()
        result[["matchtable"]] = match.table
        result[["nr.matches"]] = 0
        result[["matched_reads"]] = NA
        return(result)
    }
    
    matched.reads = unique(matches[,1])         
    nr.matching = length(matched.reads)
    cat(paste("  Number of unique matching reads with eval <=", eval, "in table '", basename(match.table), "' : ", nr.matching, "\n")) 
    if(is.null(match.fn)) match.fn = paste(match.table, "matched_reads.txt", sep="_")
    write.table(matched.reads, file = match.fn, quote = FALSE, row.names = FALSE, col.names = FALSE) 
    
    result = list()
    result[["matchtable"]] = match.table        
    result[["nr.matches"]] = nr.matching        
    result[["matched_reads"]] = match.fn          
    return(result) 
}

## overLapper- author: Thomas Girke
overLapper <- function(setlist=setlist, complexity=1:length(setlist), sep="-", cleanup=FALSE, keepdups=FALSE, type) {
	if(cleanup==TRUE) {
		setlist <- sapply(setlist, function(x) gsub("([A-Z])", "\\U\\1", x, perl=T, ignore.case=T))
		setlist <- sapply(setlist, function(x) gsub("^ {1,}| {1,}$", "", x, perl=T, ignore.case=T))
	}
	if(keepdups==TRUE) {
		dupCount <- function(setlist=setlist) {
			count <- table(setlist)
			paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".")
		}
		mynames <- names(setlist)
		setlist <- lapply(setlist, function(x) dupCount(x)) 
		names(setlist) <- mynames
	}	
	setunion <- sort(unique(unlist(setlist)))
	setmatrix <- sapply(names(setlist), function(x) setunion %in% unique(setlist[[x]])) 
	rownames(setmatrix) <- setunion
	storage.mode(setmatrix) <- "numeric"
	labels <- names(setlist)
	allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
	allcombl <- unlist(allcombl, recursive=FALSE)
	complevels <- sapply(allcombl, length)
	
	if(type=="intersects") {
		OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
		names(OLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Intersect_List=OLlist))
	}	
	if(type=="vennsets") {
		vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		}
		vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
		names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Venn_List=vennOLlist))
	}
}

vennPlot <- function(counts=counts, mymain="Venn Diagram", mysub="default", setlabels="default", yoffset=seq(0,10,by=0.34), ccol=rep(1,31), colmode=1, lcol=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), lines=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), mylwd=3, diacol=1, type="ellipse", ccex=1.0, lcex=1.0, ...) {
	if(is.list(counts)==FALSE) {
		counts <- list(counts)
	}
	
	if(!length(counts[[1]]) %in%  c(3,7,15,31)) stop("Only the counts from 2-5 way venn comparisons are supported.")
	if(length(counts[[1]])==3) {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:2], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), sep="")
		} else { 
			mysub <- mysub 
		}
		
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.1, 7.0, 5.0), y=c(6.0, 6.0, 6.0),  counts=counts[[i]])
                        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                        if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) } 
		}

		if(length(setlabels)==1 & setlabels[1]=="default") { 
			setlabels <- names(counts[[1]][1:2])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0), c(8.8, 8.8), labels=setlabels, col=lcol, cex=lcex, ...)	
	}

	if(length(counts[[1]])==7) { 
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:3], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), sep="")
		} else { 
			mysub <- mysub
		}

		symbols(x=c(4, 6, 5), y=c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.0, 7.0, 5.0, 5.0, 3.8, 6.3, 5.0), 
                                           y=c(6.5, 6.5, 3.0, 7.0, 4.6, 4.6, 5.3), 
                                           counts=counts[[i]])
	        	 if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                         if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }

		}

		if(length(setlabels)==1 & setlabels[1]=="default") { 
			setlabels <- names(counts[[1]][1:3])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels=setlabels, col=lcol, cex=lcex, ...)	
	}
	
	if(length(counts[[1]])==15 & type=="ellipse") {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:4], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
		} else { 
			mysub <- mysub
		}

		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments  
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}

		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(1.5, 3.5, 6.5, 8.5, 2.9, 3.1, 5.0, 5.0, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 5.0), 
                                                   y=c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 6.0, 2.2, 5.9, 4.0, 1.4, 1.4, 4.0, 2.8), 
                                                   counts=counts[[i]])
	        	        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                                if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE) 
		}
		ellipseVenn(...)
	} 

	if(length(counts[[1]])==15 & type=="circle") {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:4], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
		} else { 
			mysub <- mysub
		}

		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)

		for(i in seq(along=counts)) {
		        olDF <- data.frame(x=c(3.0, 6.5, 3.0, 6.5, 4.8, 3.0, 4.8, 4.8, 6.5, 4.8, 3.9, 5.7, 3.9, 5.7, 4.8), 
                                           y=c(7.2, 7.2, 3.2, 3.2, 7.2, 5.2, 0.4, 0.4, 5.2, 3.2, 6.3, 6.3, 4.2, 4.2, 5.2), 
                                           counts=counts[[i]])
			 if(colmode==1) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol, cex=ccex, ...) } # rows 14-15 of olDF are printed in next step
			 if(colmode==2) { text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol[[i]], cex=ccex[i], ...) }
			text(c(4.8), c(0.8) + yoffset[i], paste("Only in ", names(counts[[1]][1]), " & ", names(counts[[1]][4]), ": ", olDF$counts[7], "; Only in ", names(counts[[1]][2]), " & ", names(counts[[1]][3]), ": ", olDF$counts[8], sep=""), col=diacol, cex=ccex, ...)
		}

			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
		text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels=setlabels, col=lcol, cex=lcex, ...)
	} 

	if(length(counts[[1]])==31) {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:5], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), paste("; S5 =", sample_counts[5]), sep="")
		} else { 
			mysub <- mysub
		}
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments  
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			screen(1, new=FALSE)
			plotellipse(center=c(4.83,6.2), radius=c(1.43,4.11), rotate=0, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.25,5.4), radius=c(1.7,3.6), rotate=66, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.1,3.5), radius=c(1.55,3.9), rotate=150, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.48,3.15), radius=c(1.55,3.92), rotate=210, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(3.7,4.8), radius=c(1.7,3.6), rotate=293.5, segments=360, xlab="", ylab="", col=lines[5], axes=FALSE, lwd=mylwd, ...)

			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(4.85, 8.0, 7.1, 3.5, 2.0, 5.90, 4.4, 4.60, 3.60, 7.1, 6.5, 3.2, 5.4, 6.65, 3.40, 5.00, 6.02, 3.60, 5.20, 4.03, 4.20, 6.45, 6.8, 3.39, 6.03, 5.74, 4.15, 3.95, 5.2, 6.40, 5.1), 
                                                   y=c(8.30, 6.2, 1.9, 1.6, 5.4, 6.85, 6.6, 2.45, 6.40, 4.3, 6.0, 4.6, 2.1, 3.40, 3.25, 6.43, 6.38, 5.10, 2.49, 6.25, 3.08, 5.30, 4.0, 3.80, 3.20, 5.95, 5.75, 3.75, 3.0, 4.50, 4.6),
					counts=counts[[i]]) 
	        	        if(colmode==1) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...) }
                                if(colmode==2) { text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol[[i]], cex=ccex[i], ...) }
			}

			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:5])
			} else {
				setlabels <- setlabels
			}
			text(c(5.7, 7.9, 8.5, 4.2, 0.8), c(9.9, 7.9, 1.9, 0.0, 7.3), adj=c(0, 0.5), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE) 
		}
		ellipseVenn(...)
	} 
}

olBarplot <- function(OLlist=OLlist, mycol="default", margins=c(6, 10, 3, 2), mincount=0, mysub="default", ...) {
	counts <- sapply(OLlist[[4]], length)
	mylogical <- counts >= mincount
	counts <- counts[mylogical]
	
	if(mycol=="default") {
		mycol <- OLlist$Complexity_Levels
		mycol <- mycol[mylogical] 
	} else {
		mycol <- mycol	
	}

	if(mysub=="default") {
		mysub <- paste("Min Count:", mincount)
	} else {
		mysub <- mysub
	}

	par(mar=margins) # Define margins to allow long labels
	barplot(counts, col=mycol, sub=mysub, ...)
	par(mar=c(5, 4, 4, 2) + 0.1) # Set margins back to default
}  # overLapper

option_list = list(
    make_option(c("-s", "--subtrees"), type="integer", default=NULL,
            help="number of subtrees [default= %default] (MANDATORY)"),
    make_option(c("-t", "--trials"), type="integer", default=NULL,
            help="number of subsampling trials [default= %default] (MANDATORY)"),
    make_option(c("-d", "--outdir"), type="character", default=NULL,
            help="outdir [default= %default] (MANDATORY)"),
    make_option(c("-u", "--subgroup_percent"), type="character", default=NULL,
            help="outdir [default= %default] (MANDATORY)"),
    make_option(c("-o", "--venn_output"), type="character", default=NULL,
            help="outdir [default= %default] (MANDATORY)"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Define number of HMM (subtrees) and nr.trials.random.picking

nrHMM = as.integer(opt$subtrees)
nr.trials.random.picking = as.integer(opt$trials)
outfolder = opt$outdir
subgroup.percent = as.integer(opt$subgroup_percent)
venn.output = opt$venn_output 
nhmmer.evalue = 1
nrHMM
## +++ Parse nhhmer output table (the match tables obtained above)  

cat("\n  Finding number of reads matching to metagenomics probe\n\n")
hit.matrix = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking) 
rownames(hit.matrix) = paste("subtree", 1:nrHMM, sep="_") 
colnames(hit.matrix) = paste("trial", 1:nr.trials.random.picking, sep="_")

match.reads.fn = matrix(NA, nrow = nrHMM, ncol = nr.trials.random.picking)        
rownames(match.reads.fn) = paste("subtree", 1:nrHMM, sep="_")
colnames(match.reads.fn) = paste("trial", 1:nr.trials.random.picking, sep="_")

#web.df = as.data.frame(hit.matrix)        
for (i in 1:nrHMM) {
    if(0) {
    # skip this condition for the moment
    #if(!thinned[i]) {     
        match.tab = matchtbl.fn[i,1]
        match.tab = file.path(getwd(), outfolder, match.tab)
        match.reads.fn[i,1] = paste("subtree", i, "matched_reads.txt", sep="_")        
        fn = file.path(getwd(), outfolder, match.reads.fn[i,1])
        hits <- parse_nhmmer_matchtable(match.tab, match.fn = fn, eval = nhmmer.evalue)            
        hit.matrix[i,1] = hits$nr.matches
        # if(hits$nr.matches > 0) {
        #     web.df[i,1] = return.link.tag(hits$matched_reads, caption = hits$nr.matches)
        # } else {
        #     web.df[i,1] = hits$nr.matches
        # }
    } else {        
        for (k in 1:nr.trials.random.picking) {    
            # match.tab = matchtbl.fn[i,k]
            # match.tab = subtree_{n}_trial_{t}_subgroup_{subgroup_percent}_matches.tbl
            match.tab = paste("subtree", i, "trial", k, "subgroup", subgroup.percent, "matches.tbl", sep="_")
            match.tab = file.path(getwd(), outfolder, match.tab)
            match.reads.fn[i,k] = paste("subtree", i, "trial", k, "matched_reads.txt", sep="_")        
            fn = file.path(getwd(), outfolder, match.reads.fn[i,k])
            hits <- parse_nhmmer_matchtable(match.tab, match.fn = fn, eval = nhmmer.evalue)
            hit.matrix[i,k] = hits$nr.matches
            # if(hits$nr.matches > 0) {
            #     web.df[i,k] = return.link.tag(hits$matched_reads, caption = hits$nr.matches)
            # } else {
            #     web.df[i,k] = hits$nr.matches
            # }                  
        }         
    }         
}         

if(max(hit.matrix, na.rm = TRUE) == 0) {
    cat("\n\n  Sorry, not a single hit found.\n")
    #insert.line(outhtml, percent=100, size=1)
    #stop.html(outhtml)          
    cat(paste("  Intermediate results written to '", outhtml, "'\n\n"))  
}


averages = apply(hit.matrix, 1, mean, na.rm = TRUE) 
stddev = apply(hit.matrix, 1, sd, na.rm = TRUE)
coeff.var = round(stddev / averages, 2)
averages = round(averages, 2)
stddev = round(stddev, 2)  
#web.df = cbind(web.df, mean = averages, sd = stddev, cv = coeff.var) 

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Number of reads matching to the metagenome:", color="black")
# insert.table(outhtml, web.df, title="", border=1, align="center")

list.names = 1:nrHMM
matches.subtree <- vector("list", length(list.names))  
names(matches.subtree) <- list.names

list.names = 1:nr.trials.random.picking
reads.group <- vector("list", length(list.names))  
names(reads.group) <- list.names
    
for (i in 1:nrHMM) {  
    if(0) {
    # skip this condition for the moment
    #if(!thinned[i]) {     
        if(hit.matrix[i,1] > 0) {         
            fn = file.path(getwd(), outfolder, match.reads.fn[i,1])              
            if(!file.exists(fn)) stop(paste("  Final output - file '", fn, "' not found.\n\n"))
            matches.subtree[[i]] = read.table(fn, stringsAsFactors = FALSE)$V1        
        } else {
            matches.subtree[[i]] = NA 
        }                
    } else {          
        for (k in 1:nr.trials.random.picking) {    
            if(hit.matrix[i,k] > 0) {     
                fn = file.path(getwd(), outfolder, match.reads.fn[i,k])              
                if(!file.exists(fn)) stop(paste("  Final output - file '", fn, "' not found.\n\n"))
                reads.group[[k]] = read.table(fn, stringsAsFactors = FALSE)$V1  
            } else {
                reads.group[[k]] = NA  
            }                  
        }     
        for (k in 1:nr.trials.random.picking) if(!is.na(reads.group[[k]][1])) matches.subtree[[i]] = union(matches.subtree[[i]], reads.group[[k]])
    }         
}             


cat(paste("\n  Collecting the hits for the", nrHMM, "subtrees (union of all hits)\n\n")) 

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text=paste("  Collecting the hits for the", nrHMM, "subtrees (union of all hits):"), color="black")
# line.feed(outhtml, 1)
# line.feed(outhtml, 1)

for (i in 1:nrHMM) {
    if(is.na(matches.subtree[[i]][1])) {
        cat(paste("  No matching reads for subtree", i, "\n"))
        #insert.blanks(outhtml,2) 
        #insert.text(outhtml, text=paste("  Matching reads for subtree", i, ": none"), color="black")
        #line.feed(outhtml, 1)
    } else {
        fn = paste("subtree", i, "final_matches.txt", sep="_")        
        fn = file.path(getwd(), outfolder, fn)  
        write.table(matches.subtree[[i]], file = fn, quote = FALSE, row.names = FALSE, col.names = FALSE)
        cat(paste("  Matching reads for subtree", i, "saved to '", basename(fn), "'\n")) 
        # link the file:
        # insert.blanks(outhtml,2) 
        # insert.text(outhtml, text=paste("  Matching reads for subtree", i, ":"), color="black")
        # set.link(outhtml, fn, caption="txt")
        # insert.blanks(outhtml, 1)
        # insert.text(outhtml, text=paste("(", return.stained(length(matches.subtree[[i]]), "blue"), "reads )"), color="black")
        # line.feed(outhtml, 1)
    }
}

## +++ Venn diagram:

# insert.line(outhtml, percent=100, size=1)
# insert.text(outhtml, text="Comparison of the results for the subtrees:", color="black")
# line.feed(outhtml, 1)

for (i in 1:nrHMM) {
	if(is.na(matches.subtree[[i]][1])) matches.subtree[[i]] <- NULL
}

if(length(matches.subtree) < 2) {		   
	cat(paste("\n  Cannot create Venn diagram. Only", length(matches.subtree), "subtree with matches left.\n\n")) 
	warn.txt = paste("Cannot create Venn diagram. Only", length(matches.subtree), "subtree with matches left.")
	# insert.blanks(outhtml,2) 
	# put.warning(outhtml, warn.txt) 
} else {
	OLlist = overLapper(setlist = matches.subtree, sep = "_", type = "vennsets")
	counts = sapply(OLlist$Venn_List, length)
	#X11(width=6, height=6)
	#vennPlot(counts = counts)
    fn = venn.output
	outpng = fn
    outpdf = gsub("png", "pdf", fn)
    #outpng
    #fn = paste(project.id, "Venn_diagram.png", sep="_")	 
	#outpng = file.path(getwd(), outfolder, fn)
	#png(outpng, unit="px", width=960, height=960, res=120)
    pdf(outpdf)
    vennPlot(counts = counts)#, file = outpng, unit="px", width=960, height=960, res=120
    dev.off()
    #invisible(dev.copy(png, file = outpng, unit="px", width=960, height=960, res=120)); invisible(dev.off()) 
	cat(paste("\n  Venn diagram saved to '", fn, "'\n"))
	# insert.image(outhtml, outpng, width=720, border=0)
	# line.feed(outhtml, 1)
}

