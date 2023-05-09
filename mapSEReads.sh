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

#### usage ####
usage() {
	echo Program: "mapSEReads.sh (map single end reads)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: mapSEReads.sh -i <file> [OPTIONS]"
	echo "Options:"
	echo " -i <file>   [input fastq file(s) with single end reads]"
    echo "             [if multiple separate them by a comma]"
	echo "[OPTIONS]"
	echo " -m <dir>    [output directory to store mapped reads (default: .)]"
	echo " -g <string> [genome (default: mm9)]"
    echo "             [mm9, mm10, hg19, hg38, dm6, ERCC, hg19_dm6, hg19_mm9, mm9_hg19, hg38_mm10, mm10_hg38, ce11_dm6, ce11_mm10]"
    echo "             [**NOTE**: cases like mm9_hg19: assembly post '_' is considered spikeIn (hg19)]"
    echo " -p <int>    [number of processors (default: 1)]"
    echo "             [for STAR, keep it max to 20]"
    echo " -d <string> [identifier for output BAM file (default: same as fastq file)]"
    echo " -a <dir>    [directory to keep adapter trimmed fastq files]"
    echo "             [required to auto-trim adapter sequences]"
    echo " -y <dir>    [copy input fastq file(s) to specified directory for mapping]"
    echo "             [useful when working on mounted directory]"
    echo "[OPTIONS: bowtie2 (ChIP- or RNA-seq) (default)]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -U          [remove PCR duplicate reads from output bam file (samtools markdup)]"
    echo " -c          [scale the read coverage to TPM in output bigWig files (overridden if -z)]"
    echo " -C          [scale the read coverage to 1x in output bigWig files (overridden if -z)]"
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
    echo " -z <string> [spike-in normalization genome (if combined genome not available, eg. hg19_dm6)]"
    echo "             [required to auto-determine scaling factor using splike-in data]"
    echo " -Z          [only bam to bw coversion; bam file exists]"
    echo " -Y          [no bam to bw coversion]"
    echo " -r          [map reads for repeat analysis (-u option is ignored)]"
    echo " -H          [map reads for HiC data analysis (--reorder --local)]"
    echo "[OPTIONS: bwa]"
    echo " -B          [perform alignment using bwa]"
    echo "[OPTIONS: STAR (RNA-seq)]"
    echo " -S          [perform alignment accommodating for splice junctions using STAR]"
    echo " -u          [report only uniquely mapped reads]"
    echo " -c          [scale the read coverage to TPM in output bigWig files]"
    echo " -f <int>    [trim <int> bases from 5'/left end of reads (default: 0)]"
    echo " -t <int>    [trim <int> bases from 3'/right end of reads (default: 0)]"
    echo " -r          [map reads for repeat analysis (-u option is ignored)]"
    echo "[OPTIONS: Kallisto (RNA-seq)]"
    echo " -K          [perform alignment using kallisto]"
    echo " -T <int>    [average fragment length (default: 200)]"
    echo " -N <int>    [standard deviation of fragment length (default: 30)]"
    echo "[NOTE: splike-in scaling formula]"
    echo "             [https://www.sciencedirect.com/science/article/pii/S2211124714008729#app3]"
    echo "[NOTE: mapping to repeats]"
    echo "             [https://github.com/mhammell-laboratory/TEtranscripts]"
    echo "             [https://www.nature.com/scitable/topicpage/transposons-the-jumping-genes-518/]"
    echo "             [https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-019-0192-1]"
    echo "             [http://genesdev.cshlp.org/content/28/13/1410.full]"
    echo "             [https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3566-0]"
    echo "             [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1375157]"
    echo "[NOTE: mapping for HiC data]"
    echo "             [https://hicexplorer.readthedocs.io/en/latest/content/mES-HiC_analysis.html"
	echo " -h          [help]"
	echo
	exit 0
}

#### parse options ####
while getopts i:m:g:p:d:a:y:uUcCek:q:lf:t:L:I:D:E:z:ZYrHBSKT:N:h ARG; do
	case "$ARG" in
		i) FASTQ=$OPTARG;;
		m) MAPDIR=$OPTARG;;
		g) GENOME=$OPTARG;;
        p) PROCESSORS=$OPTARG;;
        d) ID=$OPTARG;;
        a) TRIM_FASTQ_DIR=$OPTARG;;
        y) COPYDIR=$OPTARG;;
        u) UNIQUE=1;;
        U) REMOVE_DUPLICATE=1;;
        c) SCALE=1;;
        C) SCALE1x=1;;
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
        z) SPIKEIN=$OPTARG;;
        Z) BAMTOBW=1;;
        Y) NOBAMTOBW=1;;
        r) REPEATS=1;;
        H) HIC=1;;
        B) BWA=1;;
        S) STAR=1;;
        K) KALLISTO=1;;
        T) KALLISTO_FL=$OPTARG;;
        N) KALLISTO_SD=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$FASTQ" -o "$HELP" ]; then
	usage
fi

## create appropriate directory structure
echo -n "Create appropriate directory structure... "
if [ ! -d "$MAPDIR" ]; then
    mkdir $MAPDIR
fi
if [ ! -z "$COPYDIR" ]; then
    echo -n "Copying fastq file(s) to $COPYDIR... "
    mkdir -p $COPYDIR
    oIFS=$IFS
    IFS=","
    FASTQFILES=($FASTQ)
    FASTQFILES_COUNT=${#FASTQFILES[@]}
    IFS=$oIFS

    FASTQ=""
    for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
        TEMPID=$(echo ${FASTQFILES[$i]} | perl -ane '$_=~s/^.*\///g; print $_;')
        if [ ! -s "$COPYDIR/$TEMPID" ]; then
            scp ${FASTQFILES[$i]} $COPYDIR/$TEMPID
        fi
        FASTQ="$FASTQ,$COPYDIR/$TEMPID"
    done
    FASTQ=$(echo $FASTQ | perl -ane '$_=~s/^\,//g; print $_;')
    MAPDIR=$COPYDIR
fi
oIFS=$IFS
IFS=","
FASTQFILES=($FASTQ)
FASTQFILES_COUNT=${#FASTQFILES[@]}
IFS=$oIFS
for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
    if [ ! -s "${FASTQFILES[$i]}" ]; then
        echo
        echo "fastq file, ${FASTQFILES[$i]} does not exist!!"
        usage
    fi
done
#echo $FASTQ; exit
echo done

echo -n "Populating files based on input genome, $GENOME (`date`).. "
if [ "$GENOME" == "mm9" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/Ensembl/NCBIM37/Bowtie2IndexWithAbundance/bowtie/Bowtie2IndexWithAbundance"
    elif [ ! -z "$STAR" ]; then
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
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/BWA/mm10.fa"
    elif [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/kallisto/Mus_musculus.GRCm38.cdna.all.idx"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Mus_musculus/mm10/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19" ]; then
    if [ ! -z "$REPENRICH" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/Ensembl/GRCh37/Bowtie2IndexInklAbundant/bowtie/genome_and_Abundant"
    elif [ ! -z "$STAR" ]; then
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
    elif [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/STAR/"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/kallisto/Homo_sapiens.GRCh38.cdna.all.idx"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/hg38/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19_ifn" ]; then
    GENOMEINDEX="/scratch/genomes/assemblies/Homo_sapiens/interferon_genes/interferon"
elif [ "$GENOME" == "dm6" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX=""
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/Dro_melanogaster/dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/Dro_melanogaster/dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ERCC" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ERCC/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ERCC/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ERCC/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg19_dm6" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "mm9_dm6" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm9_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm9_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/mm9_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi  
elif [ "$GENOME" == "hg19_mm9" -o "$GENOME" == "mm9_hg19" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_mm9/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_mm9/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg19_mm9/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "hg38_mm10" -o "$GENOME" == "mm10_hg38" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg38_mm10/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg38_mm10/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/hg38_mm10/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ce11_dm6" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_dm6/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_dm6/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_dm6/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ce11_mm10" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_mm10/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_mm10/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ce11_mm10/bowtie2/Bowtie2IndexWithAbundance"
    fi
elif [ "$GENOME" == "ecoli" ]; then
    if [ ! -z "$STAR" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ecoli/STAR"
    elif [ ! -z "$KALLISTO" ]; then
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ecoli/kallisto/"
    else
        GENOMEINDEX="/scratch/genomes/assemblies/spikeIn/ecoli/bowtie2/Bowtie2IndexWithAbundance"
    fi
else
    echo "Presently the program only support analysis for mm9, mm10, hg19, hg38, dm6, ERCC, hg19_dm6, hg19_mm9, mm9_hg19, hg38_mm10, mm10_hg38, ce11_dm6, ce11_mm10"
    echo
    usage
fi
echo done

## auto-trim adapter sequences, if required
if [ ! -z "$TRIM_FASTQ_DIR" -a -z "$BAMTOBW" ]; then
    ## parse input fastq files in an array
    oIFS=$IFS
    IFS=","
    FASTQFILES=($FASTQ)
    FASTQFILES_COUNT=${#FASTQFILES[@]}
    IFS=$oIFS

    FASTQ=""
    for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
        TEMPID=`echo ${FASTQFILES[$i]} | perl -an -F'/\,/' -e '$TEMPID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $TEMPID.=$_."_"; } $TEMPID=~s/\_$//g; print "$TEMPID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
        if [ ! -f "$TRIM_FASTQ_DIR/$TEMPID.clipped.fastq.gz" ]; then
            trimAdapters.sh -i ${FASTQFILES[$i]} -r $TRIM_FASTQ_DIR -a auto
        fi
        if [ ! -z "$COPYDIR" ]; then
            rm ${FASTQFILES[$i]}
        fi
        FASTQ="$FASTQ,$TRIM_FASTQ_DIR/$TEMPID.clipped.fastq.gz"
    done
    FASTQ=$(echo $FASTQ | perl -ane '$_=~s/^\,//g; print $_;')
fi

## retrieve file name
if [ -z "$ID" ]; then
    ID=`echo $FASTQ | perl -an -F'/\,/' -e '$ID=(); foreach(@F) { $_=~s/^.+\///g; $_=~s/\..+$//g; chomp($_); $ID.=$_."_"; } $ID=~s/\_$//g; print "$ID\n";' | perl -an -F'//' -e 'chomp($_); if(scalar(@F)>50) { $_=~s/\_R[0-9]+.*$//g; print "$_\n"; } else { print "$_\n"; }'`;
fi
FASTQ=$(echo $FASTQ | sed 's/\,/ /g')
READLENGTH=`zless $FASTQ | head -n 2 | tail -n 1 | perl -ane '$len=length($_)-1; print $len;'`;
#echo -e "$ID\t$READLENGTH"; exit;

## splike-in normalization
if [ ! -z "$SPIKEIN" ]; then
    ARGS=""
    if [ ! -z "$UNIQUE" ]; then
        ARGS="$ARGS -u";
    fi
    if [ ! -z "$EXTEND" ]; then
        ARGS="$ARGS -e";
    fi

    ID_SPIKEIN=$ID"_"$SPIKEIN

    if [ -z "$ID_SPIKEIN.mapStat" -a ! -z "$REMOVE_DUPLICATE" ]; then
        mapSEReads.sh -i $FASTQ -g $SPIKEIN -m $MAPDIR -p $PROCESSORS -d $ID_SPIKEIN -U $ARGS
        ID_SPIKEIN=$(echo $ID_SPIKEIN"_dupRemoved")
    else
        mapSEReads.sh -i $FASTQ -g $SPIKEIN -m $MAPDIR -p $PROCESSORS -d $ID_SPIKEIN $ARGS
    fi

    if [ -z "$UNIQUE" ]; then
        SCALE_SPIKEIN=$(zless $MAPDIR/$ID_SPIKEIN.mapStat | grep 'aligned exactly\|aligned >1' | perl -ane '$sum+=$F[0]; END { print "$sum\n"; }')
    else
        SCALE_SPIKEIN=$(zless $MAPDIR/$ID_SPIKEIN.mapStat | grep 'aligned exactly' | perl -ane '$sum+=$F[0]; END { print "$sum\n"; }')
    fi

    SCALE_SPIKEIN=$(echo $SCALE_SPIKEIN | perl -ane 'printf("%0.6f", 1000000/$_);');
fi
#echo $SCALE_SPIKEIN; exit

## read arguments
if [ ! -z "$STAR" ]; then
    ARGS=""
    if [ ! -z "$SCALE" ]; then
        ARGS="$ARGS --outWigNorm RPM";
    else
        ARGS="$ARGS --outWigNorm None";
    fi
    
    if [ ! -z "$UNIQUE" ]; then
        ARGS="$ARGS --outFilterMultimapNmax 1";
    fi
else
    ARGS=""
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
fi

#echo "$GENOMEINDEX $READDIR $ID"; exit;

## start analysis
if [ ! -z "$STAR" ]; then
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        FASTQ=$(echo $FASTQ | sed 's/ /\,/g')

        if [ -z "$REPEATS" ]; then
            echo "Command used: STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn <(gunzip -c $FASTQ) --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded $ARGS" >> $MAPDIR/$ID.mapStat

            STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn <(gunzip -c $FASTQ) --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded $ARGS
        else
            echo "Command used: STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn <(gunzip -c $FASTQ) --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded --outFilterMultimapNmax 10000 --winAnchorMultimapNmax 10000 --outSAMmultNmax 1 --alignTranscriptsPerReadNmax 15000 $ARGS" >> $MAPDIR/$ID.mapStat

            STAR --genomeDir $GENOMEINDEX  --runThreadN $PROCESSORS --readFilesIn <(gunzip -c $FASTQ) --outFileNamePrefix $MAPDIR/$ID --outSAMtype BAM SortedByCoordinate --clip3pNbases $TRIM3 --clip5pNbases $TRIM5 --outWigType bedGraph --outWigStrand Unstranded --outFilterMultimapNmax 10000 --winAnchorMultimapNmax 10000 --outSAMmultNmax 1 --alignTranscriptsPerReadNmax 15000 $ARGS
        fi

        mv $MAPDIR/$ID"Aligned.sortedByCoord.out.bam" $MAPDIR/$ID.bam
        zless $MAPDIR/$ID"Log.final.out" >> $MAPDIR/$ID.mapStat
        zless $MAPDIR/$ID"Log.progress.out" >> $MAPDIR/$ID.mapStat
        zless $MAPDIR/$ID"Log.out" > $MAPDIR/$ID.log
        zless $MAPDIR/$ID"SJ.out.tab" > $MAPDIR/$ID.SJ
        samtools index -@ $PROCESSORS $MAPDIR/$ID.bam
        rm $MAPDIR/$ID"Log.final.out" $MAPDIR/$ID"Log.progress.out" $MAPDIR/$ID"Log.out" $MAPDIR/$ID"SJ.out.tab"

        ## create bigwig files for visualization at the UCSC genome browser
        if [ ! -z "$UNIQUE" ]; then
            sortBed -i $MAPDIR/$ID"Signal.Unique.str1.out.bg" > $MAPDIR/$ID.bg
        else
            sortBed -i $MAPDIR/$ID"Signal.UniqueMultiple.str1.out.bg" > $MAPDIR/$ID.bg
        fi
        rm $MAPDIR/$ID"Signal.Unique.str1.out.bg"
        rm $MAPDIR/$ID"Signal.UniqueMultiple.str1.out.bg"

        ## populating files based on input genome
        if [ "$(echo $GENOME | perl -ane 'if($_=~/\_/) { print 1; } else { print 0; }')" -eq 1 ]; then
            GENOME=$(echo $GENOME | sed 's/\_.*//g')
        fi
        GENOME_FILE=$(initialize_genome -i $FINDNFRPATH/data/annotations/GENOME_FILE -g $GENOME)
        if [ ! -f "$GENOME_FILE" ]; then
            echo
            echo "computation create bigWig file for $GENOME"
            echo "please add the chromosome size file for $GENOME at $FINDNFRPATH/data/annotations"
            echo "also update the $FINDNFRPATH/data/annotations/GENOME_FILE"
            echo
        else
            bedGraphToBigWig $MAPDIR/$ID.bg $GENOME_FILE $MAPDIR/$ID.bw
        fi
    fi
elif [ ! -z "$KALLISTO" ]; then
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... "
        kallisto quant -i $GENOMEINDEX -o $MAPDIR/$ID -b 100 --single --bias -l $KALLISTO_FL -s $KALLISTO_SD -t $PROCESSORS $FASTQ
    fi
elif [ ! -z "$HIC" ]; then
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... "
        # https://hicexplorer.readthedocs.io/en/latest/content/mES-HiC_analysis.html
        bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U $FASTQ --reorder --local 2>>$MAPDIR/$ID.mapStat | samtools view -@ $PROCESSORS -Shb - > $MAPDIR/$ID.bam 
    fi
else
    if [ -z "$BAMTOBW" ]; then
        echo "Map for $ID... " >$MAPDIR/$ID.mapStat

        if [ ! -z "$BWA" ]; then
            ## command check
            echo "Command used: bwa mem $GENOMEINDEX $FASTQ -t $PROCESSORS" >>$MAPDIR/$ID.mapStat

            bwa mem $GENOMEINDEX $FASTQ -t $PROCESSORS 2>>/dev/null | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -n -m 1500M - | samtools fixmate -@ $PROCESSORS -m - $MAPDIR/$ID.bam
            samtools flagstat -@ $PROCESSORS $MAPDIR/$ID.bam >>$MAPDIR/$ID.mapStat
        else
            if [ -z "$REPEATS" ]; then
                ## command check
                echo "Command used: zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS" >>$MAPDIR/$ID.mapStat

                if [ ! -z "$UNIQUE" ]; then
                    zless $FASTQ | sed 's/.*Hendrich.*\@/@/g' | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | grep -v XS: | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -m 1500M - | samtools markdup -@ $PROCESSORS - $MAPDIR/$ID.bam
                else
                    zless $FASTQ | sed 's/.*Hendrich.*\@/@/g' | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -m 1500M - | samtools markdup -@ $PROCESSORS - $MAPDIR/$ID.bam
                fi
            else
                ## computationally not feasible
                #echo "Command used: zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -k 10000" >>$MAPDIR/$ID.mapStat
                #zless $FASTQ | sed 's/.*Hendrich.*\@/@/g' | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -k 10000 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort -m 1500M - | samtools markdup - $MAPDIR/$ID.bam
                ## limited use, no NH: tag information
                echo "Command used: zless $FASTQ | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS" >>$MAPDIR/$ID.mapStat
                zless $FASTQ | sed 's/.*Hendrich.*\@/@/g' | bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - $ALNMODE -5 $TRIM5 -3 $TRIM3 $ARGS 2>>$MAPDIR/$ID.mapStat | samtools view -@ $PROCESSORS -S -b - | samtools sort -@ $PROCESSORS -m 1500M - | samtools markdup -@ $PROCESSORS - $MAPDIR/$ID.bam
            fi
        fi

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
fi 

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

    #if [ -z "$REPEATS" ]; then
        SCALE_SPIKEIN=$(bam2spikeInScale -i $MAPDIR/${ID}_${GENOME_SPIKEIN}.bam)
        #SCALE_SPIKEIN=$(samtools flagstat $MAPDIR/$ID"_"$GENOME_SPIKEIN.bam | grep "mapped (" | cut -f 1 -d " " | perl -ane 'printf("%0.6f", 1000000/$_);');
    #else 
    #    SCALE_SPIKEIN=$(bam2spikeInScale -i $MAPDIR/$ID"_"$GENOME_SPIKEIN.bam -u);
    #fi

    SPIKEIN=$GENOME_SPIKEIN

    ID=$(echo $ID"_"$GENOME)
fi
#echo $SCALE_SPIKEIN; exit

## create bigwig files for visualization at the UCSC genome browser
if [ -z "$STAR" -a -z "$KALLISTO" -a -z "$HIC" -a -z "$NOBAMTOBW" ]; then
    if [ ! -z "$SPIKEIN" ]; then
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
echo "done"

if [ ! -z "$COPYDIR" ]; then
    echo -n "Delete copied fastq file(s) from $COPYDIR... "
    FASTQ=$(echo $FASTQ | sed 's/\,/ /g')
    rm $FASTQ
    ## parse input fastq files in an array
    # FASTQ=$(echo $FASTQ | sed 's/ /\,/g')
    # oIFS=$IFS
    # IFS=","
    # FASTQFILES=($FASTQ)
    # FASTQFILES_COUNT=${#FASTQFILES[@]}
    # IFS=$oIFS

    # for(( i=0; i<$FASTQFILES_COUNT; i++ )); do
    #    if [ ! -s "${FASTQFILES[$i]}" ]; then
    #        rm ${FASTQFILES[$i]}
    #    fi
    # done
    # echo "done"
fi

## OLD OPTIONS (TopHap2, repeats, RepEnrich): 
#elif [ ! -z "$TOPHAT2" ]; then
#    if [ -z "$BAMTOBW" ]; then
#        echo "Map for $ID... " >$MAPDIR/$ID.mapStat
#
#        echo "Command used: tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ" >>$MAPDIR/$ID.mapStat
#        tophat2 -p $PROCESSORS --b2-sensitive --transcriptome-index=$FASTAFILE --library-type=fr-unstranded -o $MAPDIR/$ID $GENOMEINDEX $FASTQ
#
#        ## compute mapping statistics
#        mv $MAPDIR/$ID/accepted_hits.bam $MAPDIR/$ID.bam
#        samtools index $MAPDIR/$ID.bam
#        zless $MAPDIR/$ID/align_summary.txt >> $MAPDIR/$ID.mapStat
#    fi
#
#    ## create bigwig files for viualization at the UCSC genome browser
#    bedtools bamtobed -i $MAPDIR/$ID.bam -bed12 | grep '^[1-9XY]' | awk '{print "chr"$0}' > $MAPDIR/$ID/accepted_hits_corrected.bed && bedtools genomecov -bg -i $MAPDIR/$ID/accepted_hits_corrected.bed -g $CHRSIZE -split > $MAPDIR/$ID/accepted_hits.bedGraph && bedGraphToBigWig $MAPDIR/$ID/accepted_hits.bedGraph $CHRSIZE $MAPDIR/$ID.bw && rm $MAPDIR/$ID/accepted_hits.bedGraph
#elif [ ! -z "$REPENRICH" ]; then
#    if [ -z "$BAMTOBW" ]; then
#        echo "Map for $ID... " >$MAPDIR/$ID.mapStat
#
#        echo "Command used: bowtie $GENOMEINDEX -p $PROCESSORS -t -m 1 -S --max $MAPDIR/$ID'_multimap.fastq' $FASTQ $MAPDIR/$ID'_unique.sam'" >> $MAPDIR/$ID.mapStat
#
#        bowtie $GENOMEINDEX -p $PROCESSORS -t -m 1 -S --max $MAPDIR/$ID"_multimap.fastq" $FASTQ $MAPDIR/$ID"_unique.sam" 2>>$MAPDIR/$ID.mapStat
#        samtools view -bS $MAPDIR/$ID"_unique.sam" > $MAPDIR/$ID"_unique.bam"
#        samtools sort $MAPDIR/$ID"_unique.bam" -o $MAPDIR/$ID"_unique_sorted.bam"
#        mv $MAPDIR/$ID"_unique_sorted.bam" $MAPDIR/$ID"_unique.bam"
#        samtools index $MAPDIR/$ID"_unique.bam"
#        #rm $MAPDIR/$ID"_unique.sam"
#    fi
#elif [ ! -z "$REPEATS" ]; then
#    if [ -z "$BAMTOBW" ]; then
#        ## inspired from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1375157
#        echo "Map for $ID... " >$MAPDIR/$ID.mapStat
#
#        ## command check
#        echo "Command used: zless $FASTQ | /home/pundhir/software/bowtie2-2.0.0-beta6/bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -M 10000" >>$MAPDIR/$ID.mapStat
#
#        zless $FASTQ | /home/pundhir/software/bowtie2-2.0.0-beta6/bowtie2 -p $PROCESSORS -x $GENOMEINDEX -U - -5 $TRIM5 -3 $TRIM3 $ARGS -D 15 -R 2 -N 0 -L 32 -i S,1,0.75 -M 10000 2>>$MAPDIR/$ID.mapStat | samtools view -S -b - | samtools sort -m 1500M - -o $MAPDIR/$ID.bam 
#
#        samtools index $MAPDIR/$ID.bam
#    fi
