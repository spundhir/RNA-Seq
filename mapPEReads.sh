#PBS -l nodes=1:ppn=4

## DEPENDENCIES
MAPDIR="."
PROCESSORS=1
GENOME="mm9"
ALNMODE="--sensitive"
FILE_ADAPTERS="/localhome/bric/xfd783/software/Trimmomatic/adapters/TruSeq3-PE.fa"
TRIM5=0
TRIM3=0
KALLISTO_FL=200
KALLISTO_SD=30
MIN_FRAGMENT_LEN=0
MAX_FRAGMENT_LEN=500
BARCODE_FILE="/tools/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt"
UMI_LEN=12

#### usage ####
usage() {
	echo Program: "mapPEReads.sh (map paired end reads)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: mapPEReads.sh -i <file> -j <file> [OPTIONS]"
	echo "Options:"
	echo " -i <file>   [input fastq file(s) with forward paired end reads]"
    echo "             [if multiple separate them by a comma]"
	echo " -j <file>   [input fastq file(s) with reverse paired end reads]"
    echo "             [if multiple separate them by a comma]"
	echo "[OPTIONS]"
	echo " -m <dir>    [output directory to store mapped reads (default: .)]"
	echo " -g <string> [genome (default: mm9)]"
    echo "             [mm9, mm10, hg19, hg38, dm6, ERCC, hg19_dm6, mm9_dm6, hg19_mm9, mm9_hg19, hg38_mm10, mm10_hg38, mm10_dm6, mm10_ecoli, ce11_dm6, ce11_mm10]"
    echo "             [**NOTE**: cases like mm9_hg19: assembly post '_' is considered spikeIn (hg19)]"
    echo "             [can also be path to directory containing STAR or BOWTIE index files]"
    echo " -p <int>    [number of processors (default: 1)]"
    echo "             [for STAR, keep it max to 20]"
    echo " -d <string> [identifier for output BAM file (default: same as fastq file)]"
    echo " -y <dir>    [copy input fastq file(s) to specified directory for mapping]"
    echo "             [useful when working on mounted directory]" 
    echo " -a          [trim adapters]"
    echo " -F <file>   [file containing adapter sequences (default: /localhome/bric/xfd783/software/Trimmomatic/adapters/TruSeq3-PE.fa)"]
    echo "[OPTIONS: bowtie2 (ChIP- or RNA-seq) (default)]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -U          [remove PCR duplicate reads from output bam file (samtools markdup)]"
    echo " -c          [scale the read coverage to TPM in output bigWig files]"
    echo " -C          [scale the read coverage to 1x in output bigWig files]"
    echo " -A <string> [scale the read coverage using CUTANA spike-in scales in output bigWig files]"
    echo "             [h3k4me1, h3k4me2, h3k4me3]"
    echo "             [h3k9me1, h3k9me2, h3k9me3]"
    echo "             [h3k27me1, h3k27me2, h3k27me3]"
    echo "             [h3k36me1, h3k36me2, h3k36me3]"
    echo "             [h3k20me1, h3k20me2, h3k20me3]"
    echo " -e          [extend 3' end of reads in output bigWig files]"
    echo " -k <int>    [instead of reporting best alignment, report input number of alignments per read]"
    echo " -q <string> [end-to-end: --very-fast, --fast, --sensitive, --very-sensitive (default: --sensitive)]"
    echo "             [local: --very-fast-local, --fast-local, --sensitive-local, --very-sensitive-local (default: --sesitive-local)"
    echo " -l          [local alignment; ends might be soft clipped]"
    echo " -f <int>    [trim <int> bases from 5'/left end of reads (default: 0)]"
    echo " -t <int>    [trim <int> bases from 3'/right end of reads (default: 0)]"
    echo " -L <int>    [length of seed substring; must be >3, <32 (default: 22)]"
    echo " -I <string> [interval between seed substrings w/r/t read len (default: S,1,1.15)]"
    echo " -D <int>    [give up extending after <int> failed extends in a row (default: 15)]"
    echo " -E <int>    [for reads w/ repetitive seeds, try <int> sets of seeds (default: 2)]"
    echo " -W <int>    [minumum fragment length (default: 0)]"
    echo " -X <int>    [maximum fragment length (default: 500)]"
    echo " -Q          [suppress unpaired alignments for paired reads (--no-mixed)]"
    echo " -R          [suppress discordant alignments for paired reads (--no-discordant)]"
    echo " -Y          [consider concordant when mates extend past each other (--dovetail)]"
    echo " -P          [suppress SAM records for unaligned reads (--no-unal)]"
    echo " -Z          [only bam to bw coversion; bam file exists]"
    echo " -r          [map reads for repeat analysis (-u option is ignored)]"
    echo "[OPTIONS: bwa]"
    echo " -B          [perform alignment using bwa]" 
    echo "[OPTIONS: STAR (RNA-seq)]"
    echo " -S          [perform alignment accommodating for splice junctions using STAR]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -c          [scale the read coverage to TPM in output bigWig files]"
    echo " -f <int>    [trim <int> bases from 5'/left end of reads (default: 0)]"
    echo " -t <int>    [trim <int> bases from 3'/right end of reads (default: 0)]"
    echo " -r          [map reads for repeat analysis (-u option is ignored)]"
    echo "[OPTIONS: STAR (single-cell RNA-seq)]"
    echo " -s          [perform alignment accommodating for splice junctions using STAR]"
    echo " -b          [barcode file (default: /tools/cellranger-7.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz)]"
    echo " -n          [length of the UMI (default: 12)]"
    echo "[OPTIONS: Kallisto (RNA-seq)]"
    echo " -K          [perform alignment using kallisto]"
    echo " -T <int>    [average fragment length (default: 200)]"
    echo " -N <int>    [standard deviation of fragment length (default: 30)]"
    echo "[NOTE: mapping single cell data]"
    echo "             [https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md]"
    echo "             [1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read]"
    echo "[NOTE: splike-in scaling formula]"
    echo "             [https://www.sciencedirect.com/science/article/pii/S2211124714008729#app3]"
    echo "[NOTE: mapping to repeats]"
    echo "             [https://github.com/mhammell-laboratory/TEtranscripts]"
    echo "             [https://www.nature.com/scitable/topicpage/transposons-the-jumping-genes-518/]"
    echo "             [https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1]"
    echo "             [http://genesdev.cshlp.org/content/28/13/1410.full]"
    echo "             [https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3566-0]"
    echo "             [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1375157]"
    echo "[NOTE: samtools flag meaning]"
    echo "             [http://broadinstitute.github.io/picard/explain-flags.html]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:j:m:g:p:d:y:aF:uUcCA:ek:q:lf:t:L:I:D:E:W:X:QRYPZrBSsb:n:KT:N:h ARG; do
	case "$ARG" in
		i) FASTQ_FORWARD=$OPTARG;;
		j) FASTQ_REVERSE=$OPTARG;;
		m) MAPDIR=$OPTARG;;
		g) GENOME=$OPTARG;;
        p) PROCESSORS=$OPTARG;;
        d) ID=$OPTARG;;
        y) COPYDIR=$OPTARG;;
        a) TRIM_ADAPTERS=1;;
        F) FILE_ADAPTERS=$OPTARG;;
        u) UNIQUE=1;;
        U) REMOVE_DUPLICATE=1;;
        c) SCALE=1;;
        C) SCALE1x=1;;
        A) SCALE_CUTANA=$OPTARG;;
        e) EXTEND=1;;
        k) ALNCOUNT=$OPTARG;;
        q) ALNMODE=$OPTARG;;
        l) LOCAL=1;;
        f) TRIM5=$OPTARG;;
        t) TRIM3=$OPTARG;;
        L) SEED=$OPTARG;;
        I) INTERVAL=$OPTARG;;
        D) GIVEUP=$OPTARG;;
        E) TRIES=$OPTARG;;
        W) MIN_FRAGMENT_LEN=$OPTARG;;
        X) MAX_FRAGMENT_LEN=$OPTARG;;
        Q) NO_MIXED=1;;
        R) NO_DISCORDANT=1;;
        Y) DOVETAIL=1;;
        P) NO_UNAL=1;;
        Z) BAMTOBW=1;;
        r) REPEATS=1;;
        B) BWA=1;;
        S) STAR=1;;
        s) STAR_SINGLE=1;;
        b) BARCODE_FILE=$OPTARG;;
        n) UMI_LEN=$OPTARG;;
        K) KALLISTO=1;;
        T) KALLISTO_FL=$OPTARG;;
        N) KALLISTO_SD=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$FASTQ_FORWARD" -o -z "$FASTQ_REVERSE" -o "$HELP" ]; then
	usage
fi

## create appropriate directory structure
echo -n "Create appropriate directory structure... "
if [ ! -z "$COPYDIR" ]; then
    echo -n "Copying fastq file(s) to $COPYDIR... "
    mkdir -p $COPYDIR
    oIFS=$IFS
    IFS=","
    FASTQFILES=($FASTQ_FORWARD)
    FASTQFILES_COUNT=${#FASTQFILES[@]}
    IFS=$oIFS

    FASTQ_FORWARD=""
    for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
        TEMPID=$(echo ${FASTQFILES[$i]} | perl -ane '$_=~s/^.*\///g; print $_;')
        if [ ! -s "$COPYDIR/$TEMPID" ]; then
            scp ${FASTQFILES[$i]} $COPYDIR/$TEMPID
        fi
        FASTQ_FORWARD="$FASTQ_FORWARD,$COPYDIR/$TEMPID"
    done
    FASTQ_FORWARD=$(echo $FASTQ_FORWARD | perl -ane '$_=~s/^\,//g; print $_;')

    oIFS=$IFS
    IFS=","
    FASTQFILES=($FASTQ_REVERSE)
    FASTQFILES_COUNT=${#FASTQFILES[@]}
    IFS=$oIFS

    FASTQ_REVERSE=""
    for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
        TEMPID=$(echo ${FASTQFILES[$i]} | perl -ane '$_=~s/^.*\///g; print $_;')
        if [ ! -s "$COPYDIR/$TEMPID" ]; then
            scp ${FASTQFILES[$i]} $COPYDIR/$TEMPID
        fi
        FASTQ_REVERSE="$FASTQ_REVERSE,$COPYDIR/$TEMPID"
    done
    FASTQ_REVERSE=$(echo $FASTQ_REVERSE | perl -ane '$_=~s/^\,//g; print $_;')
    MAPDIR=$COPYDIR
fi
if [ ! -d "$MAPDIR" ]; then
    mkdir $MAPDIR
fi
#echo "$FASTQ_FORWARD"; echo "$FASTQ_REVERSE"; exit
echo done

echo -n "Populating files based on input genome, $GENOME (`date`).. "
if [ "$GENOME" == "mm9" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie/Bowtie2IndexWithAbundance"
    elif [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm9/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm9/kallisto/Mus_musculus.NCBIM37.67.cdna.all.idx"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm9/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "mm10" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$BWA" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/bwa/mm10.fa"
    elif [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/kallisto/Mus_musculus.GRCm38.cdna.all.idx"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie/genome_and_Abundant"
    elif [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg19/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg19/kallisto/Homo_sapiens.GRCh37.cdna.all.idx"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg19/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg38" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$BWA" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/bwa/hg38.fa.chr"
    elif [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/kallisto/Homo_sapiens.GRCh38.cdna.all.idx"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19_ifn" ]; then
    GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/interferon_genes/interferon"
elif [ "$GENOME" == "dm6" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Dro_melanogaster/dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Dro_melanogaster/dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Dro_melanogaster/dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ERCC" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ERCC/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ERCC/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ERCC/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19_dm6" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "mm9_dm6" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm9_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm9_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm9_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "mm10_dm6" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm10_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm10_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm10_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "mm10_ecoli" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm10_ecoli/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm10_ecoli/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm10_ecoli/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19_mm9" -o "$GENOME" == "mm9_hg19" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_mm9/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_mm9/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_mm9/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg38_mm10" -o "$GENOME" == "mm10_hg38" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg38_mm10/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg38_mm10/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg38_mm10/bowtie2/Bowtie2IndexWithAbundance"
    fi  
elif [ "$GENOME" == "ce11_dm6" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ce11_mm10" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_mm10/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_mm10/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_mm10/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ecoli" ]; then
    if [ ! -z "$STAR" -o ! -z "$STAR_SINGLE" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ecoli/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ecoli/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ecoli/bowtie2/Bowtie2IndexWithAbundance"
    fi  
elif [ -d "${GENOME}" ]; then
    GENOMEINDEX=${GENOME}
else
    echo "Presently the program only support analysis for mm9, mm10, hg19, hg38, dm6, ERCC, hg19_dm6, mm9_dm6, mm10_dm6, mm10_ecoli, hg19_mm9, mm9_hg19, hg38_mm10, mm10_hg38, ce11_dm6, ce11_mm10"
echo
    usage
fi
echo done

## retrieve file name
if [ -z "$ID" ]; then
    ID=`echo $FASTQ_FORWARD | perl -an -F'/\,/' -e '$ID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $ID.=$_."_"; } $ID=~s/\_$//g; print "$ID\n";' | perl -an -F'//' -e 'chomp($_); $_=~s/\_R[0-9]+.*$//g; print "$_\n";'`;
fi

## trim adapter sequences
if [ "$TRIM_ADAPTERS" ]; then
    echo -n "Trim adapter sequences (`date`).. "
    FASTQDIR=$(dirname $FASTQ_FORWARD})
    FASTQFOWARDFILE=$(basename $FASTQ_FORWARD} | sed -E 's/\..*//g')
    FASTQREVERSEFILE=$(basename $FASTQ_REVERSE} | sed -E 's/\..*//g')
    if [ ! -s "${FASTQDIR}/${FASTQFOWARDFILE}_trimmed.fastq.gz" -o ! -s "${FASTQDIR}/${FASTQREVERSEFILE}_trimmed.fastq.gz" ]; then
        trimmomatic PE ${FASTQ_FORWARD} ${FASTQ_REVERSE} ${FASTQDIR}/${FASTQFOWARDFILE}_trimmed.fastq.gz ${FASTQDIR}/${FASTQFOWARDFILE}_trimmed_unpaired.fastq.gz ${FASTQDIR}/${FASTQREVERSEFILE}_trimmed.fastq.gz ${FASTQDIR}/${FASTQREVERSEFILE}_trimmed_unpaired.fastq.gz ILLUMINACLIP:${FILE_ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    fi
    FASTQ_FORWARD="${FASTQDIR}/${FASTQFOWARDFILE}_trimmed.fastq.gz"
    FASTQ_REVERSE="${FASTQDIR}/${FASTQREVERSEFILE}_trimmed.fastq.gz"
    echo done
fi
#echo "$FASTQ_FORWARD"; echo "$FASTQ_REVERSE"; exit

## read arguments
ARGS=""
if [ ! -z "$STAR" ]; then
    ARGS=""
    if [ ! -z "$SCALE" ]; then
        ARGS="$ARGS --outWigNorm RPM";
    else
        ARGS="$ARGS --outWigNorm None";
    fi

    if [ ! -z "$UNIQUE" ]; then
        ARGS="$ARGS --outFilterMultimapNmax 1"
    fi
else
    if [ ! -z "$ALNCOUNT" -a -z "$REPEATS" ]; then
        ARGS="$ARGS -k $ALNCOUNT";
    fi

    if [ ! -z "$LOCAL" ]; then
        ARGS="$ARGS --local";
    fi

    if [ ! -z "$SEED" -a -z "$REPEATS" ]; then
        ARGS="$ARGS -L $SEED";
    fi

    if [ ! -z "$INTERVAL" -a -z "$REPEATS" ]; then
        ARGS="$ARGS -i $INTERVAL";
    fi

    if [ ! -z "$GIVEUP" -a -z "$REPEATS" ]; then
        ARGS="$ARGS -D $GIVEUP";
    fi

    if [ ! -z "$TRIES" -a -z "$REPEATS" ]; then
        ARGS="$ARGS -R $TRIES";
    fi

    if [ ! -z "$NO_MIXED" -a -z "$REPEATS" ]; then
        ARGS="$ARGS --no-mixed";
    fi

    if [ ! -z "$NO_DISCORDANT" ]; then
        ARGS="$ARGS --no-discordant";
    fi

    if [ ! -z "$NO_UNAL" ]; then
        ARGS="$ARGS --no-unal";
    fi

    if [ ! -z "$DOVETAIL" ]; then
        ARGS="$ARGS --dovetail";
    fi
fi

## map reads
if [ -z "$KALLISTO" -a -z "$BAMTOBW" ]; then
    echo "Map for $ID... " >$MAPDIR/$ID.mapStat
fi

#echo "$GENOMEINDEX $READDIR $ID"; exit;

## start analysis
if [ ! -z "$STAR_SINGLE" ]; then
    echo "Command used: STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn $FASTQ_FORWARD $FASTQ_REVERSE --readFilesCommand zless --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --soloType CB_UMI_Simple --soloCBwhitelist $BARCODE_FILE --soloUMIlen $UMI_LEN --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR" >> $MAPDIR/$ID.mapStat

    STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn $FASTQ_FORWARD $FASTQ_REVERSE --readFilesCommand zless --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --soloType CB_UMI_Simple --soloCBwhitelist $BARCODE_FILE --soloUMIlen $UMI_LEN --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR --soloUMIdedup 1MM_CR

    mv $MAPDIR/$ID"Aligned.sortedByCoord.out.bam" $MAPDIR/$ID.bam
    zless $MAPDIR/$ID"Log.final.out" >> $MAPDIR/$ID.mapStat
    zless $MAPDIR/$ID"Log.progress.out" >> $MAPDIR/$ID.mapStat
    zless $MAPDIR/$ID"Log.out" > $MAPDIR/$ID.log
    zless $MAPDIR/$ID"SJ.out.tab" > $MAPDIR/$ID.SJ
    samtools index -@ $PROCESSORS $MAPDIR/$ID.bam
    rm $MAPDIR/$ID"Log.final.out" $MAPDIR/$ID"Log.progress.out" $MAPDIR/$ID"Log.out" $MAPDIR/$ID"SJ.out.tab"
elif [ ! -z "$STAR" ]; then
    if [ -z "$BAMTOBW" ]; then
        if [ -z "$REPEATS" ]; then
            echo "Command used: STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn $FASTQ_FORWARD $FASTQ_REVERSE --readFilesCommand zless --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded $ARGS" >> $MAPDIR/$ID.mapStat

            STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn $FASTQ_FORWARD $FASTQ_REVERSE --readFilesCommand zless --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded $ARGS
        else
            echo "Command used: STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn $FASTQ_FORWARD $FASTQ_REVERSE --readFilesCommand zless --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded --outFilterMultimapNmax 10000 --winAnchorMultimapNmax 10000 --outSAMmultNmax 1 --alignTranscriptsPerReadNmax 15000 $ARGS" >> $MAPDIR/$ID.mapStat

            STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn $FASTQ_FORWARD $FASTQ_REVERSE --readFilesCommand zless --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded --outFilterMultimapNmax 10000 --winAnchorMultimapNmax 10000 --outSAMmultNmax 1 --alignTranscriptsPerReadNmax 15000 $ARGS
        fi

        mv $MAPDIR/$ID"Aligned.sortedByCoord.out.bam" $MAPDIR/$ID.bam
        zless $MAPDIR/$ID"Log.final.out" >> $MAPDIR/$ID.mapStat
        zless $MAPDIR/$ID"Log.progress.out" >> $MAPDIR/$ID.mapStat
        zless $MAPDIR/$ID"Log.out" > $MAPDIR/$ID.log
        zless $MAPDIR/$ID"SJ.out.tab" > $MAPDIR/$ID.SJ
        samtools index -@ $PROCESSORS $MAPDIR/$ID.bam
        rm $MAPDIR/$ID"Log.final.out" $MAPDIR/$ID"Log.progress.out" $MAPDIR/$ID"Log.out" $MAPDIR/$ID"SJ.out.tab"

        ## create bigwig files for visualization at the UCSC genome browser
        # ls *.bg | parallel -j 1 'sort -k 1,1 -k 2n,2 {} > {.}.sort'
        # ls *.sort | parallel -j 1 'bedGraphToBigWig {} /scratch/genomes/annotations/SIZE/mouse.mm9.genome bigWig/{.}.bw' &
        if [ ! -z "$UNIQUE" ]; then
            sortBed -i $MAPDIR/$ID"Signal.Unique.str1.out.bg" > $MAPDIR/$ID.bg
        else
            sortBed -i $MAPDIR/$ID"Signal.UniqueMultiple.str1.out.bg" > $MAPDIR/$ID.bg
        fi
        rm $MAPDIR/$ID"Signal.Unique.str1.out.bg"
        rm $MAPDIR/$ID"Signal.UniqueMultiple.str1.out.bg"

        ## populating files based on input genome
        if [ ! -d "${GENOME}" ]; then
            if [ "$(echo $GENOME | perl -ane 'if($_=~/\_/) { print 1; } else { print 0; }')" -eq 1 ]; then
                GENOME=$(echo $GENOME | sed 's/\_.*//g')
            fi
            GENOME_FILE=$(initialize_genome -i $FINDNFRPATH/data/annotations/GENOME_FILE -g $GENOME)
            if [ ! -f "$GENOME_FILE" ]; then
                echo
                echo "Cannot create bigWig file for $GENOME"
                echo "please add the chromosome size file for $GENOME at $FINDNFRPATH/data/annotations"
                echo "also update the $FINDNFRPATH/data/annotations/GENOME_FILE"
                echo
            else
                bedGraphToBigWig $MAPDIR/$ID.bg $GENOME_FILE $MAPDIR/$ID.bw
                rm $MAPDIR/$ID.bg
            fi
        fi
    fi
elif [ ! -z "$KALLISTO" ]; then
    echo "Map for $ID... "
    kallisto quant -i $GENOMEINDEX -o $MAPDIR/$ID -b 100 --bias -l $KALLISTO_FL -s $KALLISTO_SD -t $PROCESSORS $FASTQ_FORWARD $FASTQ_REVERSE
else
    if [ -z "$BAMTOBW" ]; then
        if [ ! -z "$BWA" ]; then
            ## command check
            echo "Command used: bwa mem $GENOMEINDEX $FASTQ_FORWARD $FASTQ_REVERSE -t $PROCESSORS" >>$MAPDIR/$ID.mapStat

            bwa mem $GENOMEINDEX $FASTQ_FORWARD $FASTQ_REVERSE -t $PROCESSORS 2>>/dev/null | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -n -m 1500M - | samtools fixmate -@ $PROCESSORS -m - $MAPDIR/$ID.bam
            samtools flagstat -@ $PROCESSORS $MAPDIR/$ID.bam >>$MAPDIR/$ID.mapStat
        else
            if [ -z "$REPEATS" ]; then
                ## command check
                echo "Command used: bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS" >>$MAPDIR/$ID.mapStat

                if [ ! -z "$UNIQUE" ]; then
                    bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS 2>>$MAPDIR/$ID.mapStat | grep -v XS: | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -n -m 1500M - | samtools fixmate -@ $PROCESSORS -m - $MAPDIR/$ID.bam
                else
                    bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS 2>>$MAPDIR/$ID.mapStat | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -n -m 1500M - | samtools fixmate -@ $PROCESSORS -m - $MAPDIR/$ID.bam
                fi
     
                #alignment (options...) \
                #| samtools fixmate -m - - \
                #| samtools sort -O BAM \
                #| tee out_sorted_withDuplicates.bam \
                #| samtools markdup -r - out_sorted_Markdup.bam
            else
                ## computationally not feasible
                #echo "Command used: bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -k 10000 --no-mixed" >>$MAPDIR/$ID.mapStat
                #bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -k 10000 --no-mixed 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort -n -m 1500M - | samtools fixmate -m - $MAPDIR/$ID.bam
                ## limited use, no NH: tag information
                echo "Command used: bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS" >>$MAPDIR/$ID.mapStat
                bowtie2 -1 $FASTQ_FORWARD -2 $FASTQ_REVERSE -p $PROCESSORS -x $GENOMEINDEX $ALNMODE -5 $TRIM5 -3 $TRIM3 -I $MIN_FRAGMENT_LEN -X $MAX_FRAGMENT_LEN $ARGS 2>>$MAPDIR/$ID.mapStat | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -n -m 1500M - | samtools fixmate -@ $PROCESSORS -m - $MAPDIR/$ID.bam
            fi
        fi

        TMP=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
        samtools sort -@ $PROCESSORS $MAPDIR/$ID.bam | samtools markdup -@ $PROCESSORS - $MAPDIR/$ID"_"$TMP.bam
        mv $MAPDIR/$ID"_"$TMP.bam $MAPDIR/$ID.bam

        ## remove PCR duplicate reads (important for ChIP-seq data)
        if [ ! -z "$REMOVE_DUPLICATE" ]; then
            samtools markdup -@ $PROCESSORS -r $MAPDIR/$ID.bam $MAPDIR/$ID"_dupRemoved".bam
            ID=$(echo $ID"_dupRemoved")
        fi

        ## compute mapping statistics
        ## idxstat format: The output is TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 
        #samtools index $MAPDIR/$ID.bam && samtools idxstats $MAPDIR/$ID.bam > $MAPDIR/$ID.MappingStatistics.txt && perl -ane 'print "$F[0]\t$F[2]\t'$ID'\n";' $MAPDIR/$ID.MappingStatistics.txt >> $MAPDIR/concatenated_accepted_MappingStatistics.txt
        samtools index -@ $PROCESSORS $MAPDIR/$ID.bam
    else
        if [ ! -z "$REMOVE_DUPLICATE" ]; then
            ID=$(echo $ID"_dupRemoved")
        fi
    fi
<<"COMMENT"
COMMENT

    if [ ! -d "${GENOME}" ]; then
        ## split bam file, if mapped to multiple genome assemblies (splike-in control)
        if [ "$(echo $GENOME | perl -ane 'if($_=~/\_/) { print 1; } else { print 0; }')" -eq 1 ]; then
            GENOME_SPIKEIN=$(echo $GENOME | sed 's/.*\_//g')

            CHROM=$(samtools view -@ $PROCESSORS -H $MAPDIR/$ID.bam | awk '$1 == "@SQ" {sub("SN:", "", $2); print $2}' | grep "_"$GENOME_SPIKEIN | perl -ane 'chomp($_); print "$_ ";');
            #echo "$GENOME_SPIKEIN $CHROM"; exit

            if [ -z "$BAMTOBW" ]; then
                ## samtools view misses a lot of read in the output as happened in MLL-AF9 project
                #samtools view -S -h $MAPDIR/$ID.bam $CHROM | sed 's/_'$GENOME_SPIKEIN'//g' | samtools view -b - > $MAPDIR/$ID"_"$GENOME_SPIKEIN.bam
                cat <(samtools view -@ $PROCESSORS -H $MAPDIR/$ID.bam) <(samtools view -@ $PROCESSORS -S $MAPDIR/$ID.bam | grep "_${GENOME_SPIKEIN}") | sed 's/_'$GENOME_SPIKEIN'//g' | samtools view -@ $PROCESSORS -b - > $MAPDIR/$ID"_"$GENOME_SPIKEIN.bam
            fi

            GENOME=$(echo $GENOME | sed 's/\_.*//g')

            CHROM=$(samtools view -@ $PROCESSORS -H $MAPDIR/$ID.bam | awk '$1 == "@SQ" {sub("SN:", "", $2); print $2}' | grep "_"$GENOME | perl -ane 'chomp($_); print "$_ ";');

            if [ -z "$BAMTOBW" ]; then
                ## samtools view misses a lot of read in the output as happened in MLL-AF9 project
                #samtools view -S -h $MAPDIR/$ID.bam $CHROM | sed 's/_'$GENOME'//g' | samtools view -b - > $MAPDIR/$ID"_"$GENOME.bam
                cat <(samtools view -@ $PROCESSORS -H $MAPDIR/$ID.bam) <(samtools view -@ $PROCESSORS -S $MAPDIR/$ID.bam | grep "_${GENOME}") | sed 's/_'$GENOME'//g' | samtools view -@ $PROCESSORS -b - > $MAPDIR/$ID"_"$GENOME.bam

                samtools index -@ $PROCESSORS $MAPDIR/$ID"_"$GENOME.bam
            fi

            SCALE_SPIKEIN=$(bam2spikeInScale -i $MAPDIR/${ID}_${GENOME_SPIKEIN}.bam)
            #SCALE_SPIKEIN=$(samtools flagstat $MAPDIR/$ID"_"$GENOME_SPIKEIN.bam | grep "properly paired (" | cut -f 1 -d " " | perl -ane '$_=sprintf("%0.0f", $_/2); printf("%0.6f", 1000000/$_);');

            ID=$(echo $ID"_"$GENOME)
        elif [ ! -z "${SCALE_CUTANA}" ]; then
            SCALE_SPIKEIN=$(cutana.sh -i $FASTQ_FORWARD,$FASTQ_REVERSE -H ${SCALE_CUTANA} | cut -f 4 | perl -ane '$_=sprintf("%0.0f", $_/2); $sum+=$_; END { print "$sum"; }' | perl -ane 'printf("%0.6f", 1000000/$_);');
            if [ -z "$SCALE_SPIKEIN" ]; then
                >&2 echo -e "Error: cannot compute CATANA spike-in scale for ${SCALE_CUTANA} using $FASTQ_FORWARD,$FASTQ_REVERSE"
                >&2 echo -e "Command used: cutana.sh -i $FASTQ_FORWARD,$FASTQ_REVERSE -H ${SCALE_CUTANA}"
                exit
            fi
        fi
        #echo "$SCALE_SPIKEIN"; exit

        ## create bigwig files for visualization at the UCSC genome browser
        if [ ! -z "$SCALE_SPIKEIN" ]; then
            echo -e "Using spike-in scale, $SCALE_SPIKEIN to normalize bigWig files.. "
            if [ ! -z "$EXTEND" ]; then
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -z $SCALE_SPIKEIN -p $PROCESSORS
            else
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -z $SCALE_SPIKEIN -p $PROCESSORS
            fi
            echo "done"
        elif [ ! -z "$SCALE" ]; then
            echo -e "Using RPKM scaling to normalize bigWig files.. "
            if [ ! -z "$EXTEND" ]; then
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -c -p $PROCESSORS
            else
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -c -p $PROCESSORS
            fi
            echo "done"
        elif [ ! -z "$SCALE1x" ]; then
            echo -e "Using 1x scaling to normalize bigWig files.. "
            if [ ! -z "$EXTEND" ]; then
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -C -p $PROCESSORS
            else
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -C -p $PROCESSORS
            fi
            echo "done"
        else
            echo -e "Using no scaling to normalize bigWig files.. "
            if [ ! -z "$EXTEND" ]; then
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -e -p $PROCESSORS
            else
                bam2bwForChIP -i $MAPDIR/$ID.bam -o $MAPDIR/ -g $GENOME -p $PROCESSORS
            fi
            echo "done"
        fi
    fi
fi 

echo "done"

if [ ! -z "$COPYDIR" ]; then
    echo -n "Delete copied fastq file(s) from $COPYDIR... "
    FASTQ_FORWARD=$(echo $FASTQ_FORWARD | sed 's/\,/ /g')
    FASTQ_REVERSE=$(echo $FASTQ_REVERSE | sed 's/\,/ /g')
    rm $FASTQ_FORWARD $FASTQ_REVERSE
    echo "done"
fi
