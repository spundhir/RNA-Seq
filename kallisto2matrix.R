#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--configFile"), help="configuration file containing kallisto mapped file information (format: <id> <dir> | can be stdin)"),
	make_option(c("-o", "--outFile"), help="output file name (name.gene_counts and name.transcript_counts will be created)"),
	make_option(c("-g", "--genome"), default="mm9", help="organism genome (default: %default; Options: mm9, mm10, hg19, hg38)"),
	make_option(c("-t", "--tpm"), help="create tpm matrix (default: count)", action="store_true"),
	make_option(c("-s", "--sleuth"), help="create count matrix after batch correction using sleuth (default: count)", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$configFile) | is.null(opt$outFile)) {
	cat("\nProgram: kallisto2matrix.R (R script to compute gene count matrix from kallisto results)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))

## read configuration file
if(identical(opt$configFile, "stdin")==T) {
    configFile <- read.table(file("stdin"), header=F)
} else {
    configFile <- read.table(opt$configFile, header=F)
}
#configFile <- configFile[order(configFile$V1),]

################################################## BIOMART USAGE
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

## read each file into a data frame
abundance.lst <- apply(configFile, 1, function(x) read.table(sprintf("%s/abundance.tsv", x[2]), header=T))
abundance.df <- do.call(cbind.data.frame, abundance.lst)

##################################
## determine transcript expression
##################################
if(is.null(opt$tpm)) { 
    ## retrieve columns containing counts
    transcript_expr <- abundance.df[,c(which(names(abundance.df) %in% "est_counts"))]
} else {
    ## retrieve columns containing tpm
    transcript_expr <- abundance.df[,c(which(names(abundance.df) %in% "tpm"))]
}

transcript_expr <- cbind(abundance.df$target_id, round(transcript_expr))
colnames(transcript_expr) <- c("target_id", as.vector(configFile$V1))
transcript_expr$target_id <- gsub("\\..*", "", transcript_expr$target_id)

## add transcript coordinate information
coor <- getBM(filters="ensembl_transcript_id", attributes=c("ensembl_transcript_id", "chromosome_name", "start_position", "end_position", "strand"), values=transcript_expr$target_id, mart=mart)
coor$chromosome_name <- sprintf("chr%s", coor$chromosome_name)
coor[which(coor$strand=="1"),]$strand <- "+"
coor[which(coor$strand=="-1"),]$strand <- "-"
coor$score <- 0
coor$Length <- coor$end_position - coor$start_position
coor$Copies <- NA
result <- merge(coor, transcript_expr, by.x="ensembl_transcript_id", by.y="target_id")
result <- result[,c(c("chromosome_name", "start_position", "end_position", "ensembl_transcript_id", "score", "strand", "Length", "Copies"), colnames(result)[9:ncol(result)])]
result <- result[which(grepl("[\\_\\.]+", result$chromosome_name)==F),]
colnames(result)[colnames(result) == "ensembl_transcript_id"] <- "name"

if(is.null(opt$tpm)) {
    write.table(result, sprintf("%s.transcripts_count", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
} else {
    write.table(result, sprintf("%s.transcripts_tpm", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
}

##################################
## determine gene expression
##################################
## determine external gene name for ensembl id
tx2gene <- getBM(filters="ensembl_transcript_id", attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), values=transcript_expr$target_id, mart=mart)
files <- apply(configFile, 1, function(x) sprintf("%s/abundance.tsv", x[2]))

if(is.null(opt$tpm)) { 
    ## collapse transcript counts into gene counts
    gene_expr <- tximport(files, type="kallisto", tx2gene=tx2gene[,c(1,3)], ignoreTxVersion=T)
} else {
    ## collapse transcript tpm into gene tpm
    gene_expr <- tximport(files, type="kallisto", tx2gene=tx2gene[,c(1,3)], countsFromAbundance="lengthScaledTPM", ignoreTxVersion=T)
}

gene_expr <- cbind(as.data.frame(rownames(gene_expr$counts)), round(gene_expr$counts))
colnames(gene_expr) <- c("target_id", as.vector(configFile$V1))

## add gene coordinate information
coor <- getBM(filters="external_gene_name", attributes=c("external_gene_name", "chromosome_name", "start_position", "end_position", "strand"), values=gene_expr$target_id, mart=mart)
coor$external_gene_name <- gsub("\\..*", "", coor$external_gene_name)
coor$chromosome_name <- sprintf("chr%s", coor$chromosome_name)
coor[which(coor$strand=="1"),]$strand <- "+"
coor[which(coor$strand=="-1"),]$strand <- "-"
coor$score <- 0
coor$Length <- coor$end_position - coor$start_position
coor$Copies <- NA
result <- merge(coor, gene_expr, by.x="external_gene_name", by.y="target_id")
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

if(is.null(opt$tpm)) { 
    write.table(result, sprintf("%s.genes_count", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
} else {
    write.table(gene_tpm, sprintf("%s.genes_tpm", opt$outFile), sep="\t", quote=F, col.names=T, row.names=F)
}

q()
