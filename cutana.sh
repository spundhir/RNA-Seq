#PBS -l nodes=1:ppn=4

## DEPENDENCIES
MAPDIR="."
PROCESSORS=1
GENOME="mm9"
ALNMODE="--sensitive"
TRIM5=0
TRIM3=0
KALLISTO_FL=200
KALLISTO_SD=30

BARCODE_FILE="/localhome/bric/xfd783/software/myScripts/RNAseq/BARCODES_CUTANA"

#### usage ####
usage() {
	echo Program: "cutana.sh (determine scaling factor for CUTANA ChIP-seq protocol)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: cutana.sh -i <file> [OPTIONS]"
	echo "Options:"
	echo " -i <file>   [input fastq file(s)]"
    echo "             [if multiple (R1 & R2) separate them by a comma]"
	echo "[OPTIONS]"
	echo " -b <file>   [barcode file (format: barcode id; default: ~/software/myScripts/RNAseq/BARCODES_CUTANA)]"
    echo " -H <string> [show count for the given histone modification (default: all)]"
    echo "             [unmodified]"
    echo "             [h3k4me1, h3k4me2, h3k4me3]"
    echo "             [h3k9me1, h3k9me2, h3k9me3]"
    echo "             [h3k27me1, h3k27me2, h3k27me3]"
    echo "             [h3k36me1, h3k36me2, h3k36me3]"
    echo "             [h3k20me1, h3k20me2, h3k20me3]"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:b:H:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
        b) BARCODE_FILE=$OPTARG;;
        H) HISTONE_MOD=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$FASTQ" -o "$HELP" ]; then
	usage
fi

## organize fastq files and extract names
NAME=${FASTQ}
FASTQ=$(echo $FASTQ | sed -E 's/\,/ /g')
if [ "$(for i in $FASTQ; do echo ${i} | sed -E 's/\_R.*//g'; done | sort | uniq | wc -l)" -eq 1 ]; then
    NAME=$(for i in $FASTQ; do echo ${i} | sed -E 's/\_R.*//g'; done | sort | uniq)
fi

## start analysis
>&2 echo -e "Counting number of CUTANA barcodes for ${NAME} (`date`).. "
if [ ! -z "${HISTONE_MOD}" ]; then
    grep -v "#" $BARCODE_FILE | grep -vE "^\s+$" | grep -w ${HISTONE_MOD} | while read -a line; do
        N=$(zgrep -c ${line[0]} ${FASTQ} | cut -f 2 -d ":" | awk '{s+=$1} END {print s}');
        echo -e "${NAME}\t${line[1]}\t${line[0]}\t${N}"
    done
else
    grep -v "#" $BARCODE_FILE | grep -vE "^\s+$" | while read -a line; do
        N=$(zgrep -c ${line[0]} ${FASTQ} | cut -f 2 -d ":" | awk '{s+=$1} END {print s}');
        echo -e "${NAME}\t${line[1]}\t\t${line[0]}\t${N}"
    done
fi
>&2 echo "done"
