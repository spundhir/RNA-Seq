#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("seqLogo"))
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("session"))

## parse command line arguments
option_list <- list(
    make_option(c("-i", "--input"), help="input Position Weight Matrix (can be stdin)"),
    make_option(c("-o", "--output"), help="output seqlogo pdf file (optional)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$input))) {
    cat("\nProgram: pwm2seqlogo.R (plot sequence logo corresponding to input PWM\n")
    cat("Author: BRIC, University of Copenhagen, Denmark\n")
    cat("Version: 1.0\n")
    cat("Contact: pundhir@binf.ku.dk\n");
    print_help(parser)
    q()
}

## credits for the code goes to https://gist.github.com/dianalow/9c9a3b1beed3367300d5
if(identical(opt$input, "stdin")==T) {
    pwm <- read.table(file("stdin"), row.names = NULL, colClasses="numeric")
} else {
    pwm <- read.table(opt$input, row.names=NULL, colClasses="numeric")
}
colnames(pwm) <- c("A", "C", "G", "T")

#plot seq logo
plot_seqlogo <-function(pwm_arg){
    require(seqLogo)
    seqLogo(t(pwm_arg),ic.scale=FALSE) #ic.scale determines either frequency or bits
    #seqLogo(makePWM(pwm),ic.scale=TRUE) #ic.scale determines either frequency or bits
}

if(is.null(opt$output)) {
    outfile=sprintf("%s.pdf", opt$input)
} else {
    outfile=opt$output
}

pdf(outfile, width=40)
plot_seqlogo(pwm)
not_required <- dev.off()

save.session("test.session")
