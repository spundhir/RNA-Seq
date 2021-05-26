#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--configFile"), help="input configuration file containing bam file information (can be stdin; FORMAT: tpm test.bam)"),
	make_option(c("-o", "--outFile"), help="output file containing read counts"),
    make_option(c("-j", "--gtfFile"), help="input file containing genomic coordinate of genes in GTF format (optional)"),
    make_option(c("-n", "--genome"), default="mm9", help="genome (default: %default; mm9, mm10, hg19, hg38)"),
    make_option(c("-t", "--featureType"), default="exon", help="specify the feature type (default: %default)"),
    make_option(c("-r", "--attributeType"), default="gene_name", help="specify the feature type (default: %default)"),
    make_option(c("-O", "--allowMultiassign"), action="store_true", help="allow reads assignment to multiple features"),
    make_option(c("-M", "--countMultimapping"), action="store_true", help="also count multi-mapping reads/fragments"),
    make_option(c("-s", "--strandSpecific"), default=0, help="0: unstranded, 1: stranded, 2: reversely stranded (default: %default)"),
    make_option(c("-T", "--processors"), default=1, help="number of processors (default: %default)"),
    make_option(c("-F", "--fraction"), action="store_true", help="If specified, a fractional count 1/n will be generated for each multi-mapping read"),
    make_option(c("-p", "--primaryOnly"), action="store_true", help="Count primary alignments only. Primary alignments are identified using bit 0x100 in SAM/BAM FLAG field"),
    make_option(c("-d", "--ignoreDup"), action="store_true", help="Ignore duplicate reads in read counting. Duplicate reads are identified using bit Ox400 in BAM/SAM FLAG field"),
    make_option(c("-J", "--juncCounts"), action="store_true", help="Count number of reads supporting each exon-exon junction. Junctions were identified from those exon-spanning reads in the input (containing 'N' in CIGAR string)"),
    make_option(c("-P", "--isPairedEnd"), action="store_true", help="If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads"),
    make_option(c("-R", "--repeats"), action="store_true", help="analyze for repeats (no genomic coordinate information)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$configFile) | is.null(opt$outFile) | (is.null(opt$gtfFile) & is.null(opt$genome))) {
	cat("\nProgram: gtf2expr.R (R script to compute read count corresponding to input GTF file)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library("Rsubread"))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(session))

if(identical(opt$configFile, "stdin")==T) { 
    data <- read.table(file("stdin"), stringsAsFactors=F)
} else {
    data <- read.table(opt$configFile, stringsAsFactors=F)
}

################################################
## BIOMART USAGE
# listMarts(); ensembl = useMart("ENSEMBL_MART_ENSEMBL"); listDatasets(ensembl); listAttributes(mart)
# Use the http://www.ensembl.org website and go down the bottom of the page. Click on ’view in Archive’ and select the archive you need. Copy the url and use that url as shown below to connect to the specified BioMart database. The example below shows how to query Ensembl 54.
# > listMarts(host = "may2009.archive.ensembl.org")
# > ensembl54 = useMart(host = "may2009.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
################################################

## NOTE: external_gene_name is only available with host=grch37.ensembl.org
# listEnsembl(); listMarts(); listDatasets()
if(opt$genome=="hg19") {
    mart = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$genome=="mm9") {
    mart = useMart(host = "may2012.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
} else if(opt$genome=="hg38") {
    mart = useMart(host = "apr2020.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$genome=="mm10") {
    mart = useMart(host = "apr2020.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
} else {
    cat("Unknown genome provided\n");
    print_help(parser)
    q() 
}

## determine gene expression counts
if(is.null(opt$gtfFile)) {
    counts <- featureCounts(data$V2, annot.inbuilt=opt$genome, isGTFAnnotationFile=T, GTF.featureType=opt$featureType, GTF.attrType=opt$attributeType, allowMultiOverlap=!is.null(opt$allowMultiassign), countMultiMappingReads=!is.null(opt$countMultimapping), strandSpecific=opt$strandSpecific, nthreads=opt$processors, fraction=!is.null(opt$fraction), primaryOnly=!is.null(opt$primaryOnly), ignoreDup=!is.null(opt$ignoreDup), juncCounts=!is.null(opt$juncCounts), isPairedEnd=!is.null(opt$isPairedEnd))
} else {
    counts <- featureCounts(data$V2, annot.ext=opt$gtfFile, isGTFAnnotationFile=T, GTF.featureType=opt$featureType, GTF.attrType=opt$attributeType, allowMultiOverlap=!is.null(opt$allowMultiassign), countMultiMappingReads=!is.null(opt$countMultimapping), strandSpecific=opt$strandSpecific, nthreads=opt$processors, fraction=!is.null(opt$fraction), primaryOnly=!is.null(opt$primaryOnly), ignoreDup=!is.null(opt$ignoreDup), juncCounts=!is.null(opt$juncCounts), isPairedEnd=!is.null(opt$isPairedEnd))
}

## many gene names do not have unique genomic coordinates (for example, 5S_rRNA)
## meaning that they have multiple copies spread across different chromosomes
## table(unlist(lapply(all$Chr, function(x) length(unique(unlist(strsplit(x, ";")))))))
counts$annotation$Copies <- unlist(lapply(counts$annotation$Chr, function(x) length(unique(unlist(strsplit(x, ";"))))))

## normalize read counts by spikeIn scale factors
if(is.numeric(data$V3)==T) {
    counts$spikeInScaled <- (t(t(counts$counts) * data$V3))
    all <- cbind(counts$annotation[,c("GeneID", "Length", "Copies")], round(counts$spikeInScaled))
} else {
    all <- cbind(counts$annotation[,c("GeneID", "Length", "Copies")], round(counts$counts))
}

if(is.null(opt$repeats)) {
    ## add gene coordinate information
    coor <- getBM(filters="external_gene_name", attributes=c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand"), values=all$GeneID, mart=mart)
    #coor$external_gene_name <- gsub("\\..*", "", coor$external_gene_name)
    coor$chromosome_name <- sprintf("chr%s", coor$chromosome_name)
    coor[which(coor$strand=="1"),]$strand <- "+"
    coor[which(coor$strand=="-1"),]$strand <- "-"
    coor$score <- 0
    result <- merge(coor, all, by.x="external_gene_name", by.y="GeneID")
    result <- result[,c(c("chromosome_name", "start_position", "end_position", "external_gene_name", "score", "strand", "Length", "Copies"), colnames(result)[9:ncol(result)])]

    ## select highest expressed version of multiple copies of a gene (for example: CKS1B in hg19)
    result <- result[which(grepl("[\\_\\.]+", result$chromosome_name)==F),]
    if(!is.null(ncol(result[,c(9:ncol(result))]))) { 
        result$mean <- rowMeans(result[,c(9:ncol(result))])
        result <- result[order(result$external_gene_name, -result$mean),]
    } else {
        result <- result[order(result$external_gene_name, -result[,9]),]
    }
    result <- result[match(unique(result$external_gene_name), result$external_gene_name),]
    result$mean <- NULL
    colnames(result)[colnames(result) == "external_gene_name"] <- "name"

    write.table(result, opt$outFile, sep="\t", quote = F, col.names = T, row.names = F)
} else {
    colnames(all)[colnames(all) == "GeneID"] <- "name"

    write.table(all, opt$outFile, sep="\t", quote = F, col.names = T, row.names = F)
}
save.session("featureCounts.session")
