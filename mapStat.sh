#!/bin/bash
#PBS -l nodes=1:ppn=4

#### usage ####
usage() {
	echo Program: "mapStat (tabulate mapping statistics)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: mapStat.sh -s <file>"
	echo " -s <file>   [input mapStat file from map[S|P]EReads.sh (required)]"
	echo "[Options]"
    echo " -m <file>   [input BAM file containing mapped reads (optional)]"
    echo "             [will be used to compute additional mapping statistics]"
	echo " -q <file>   [input FASTQ (raw reads) file (optional)]"
	echo " -f <file>   [input FASTQ (clipped read) file (optional)]"
    echo " -S          [compute mapping statistics for STAR alignment results]"
	echo " -h          [help]"
    echo "[NOTE]"
    echo "             [mapping statistics are computed using:]"
    echo "             [1. .mapStat file (from bowtie2) -> id; #reads (for mapping); #reads (paired); #reads (unmapped); #reads (aligned 1 time); #reads (aligned >1 time)]"
    echo "             [2. samtools idxstats (spike-in only) -> dm6 (mapped); mm9 (mapped)]"
    echo "             [3. samtools flagstat -> #reads (QC-passed); #reads (mapped); #reads (paired); #reads (singleton); #reads (PCR duplicates)]"
	echo
	exit 0
}

MAPPING_FREQUENCY=1

#### parse options ####
while getopts s:m:q:f:Sh ARG; do
	case "$ARG" in
		s) MAPSTATFILE=$OPTARG;;
        m) BAMFILE=$OPTARG;;
		q) RAWFASTQFILE=$OPTARG;;
		f) CLIPPEDFASTQFILE=$OPTARG;;
        S) STAR=1;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ ! -f "$MAPSTATFILE" -o "$HELP" ]; then
	usage
fi

## create index of input bam file, if does not exist
if [ ! -z "$BAMFILE" -a ! -f "$BAMFILE.bai" ]; then
    samtools index $BAMFILE
fi

## bowtie2-based results
if [ -z "$STAR" ]; then
    ## compute mapping statistics for single-end data
    if [ "$(zless $MAPSTATFILE| grep -w "paired" | wc -l)" -eq 0 ]; then
        ## count total number of raw reads
        RAW_READS_COUNT="NA"
        if [ ! -z "$RAWFASTQFILE" -a -f "$RAWFASTQFILE" ]; then
            RAW_READS_COUNT=$(zless $RAWFASTQFILE | grep "@" | wc -l)
        fi

        ## count total number of reads left after quality check
        CLIPPED_READS_COUNT="NA"
        PER="NA"
        if [ ! -z "$CLIPPED_READS_COUNT" -a -f "$CLIPPEDFASTQFILE" ]; then
            CLIPPED_READS_COUNT=$(zless $CLIPPEDFASTQFILE | grep "@" | wc -l)
            if [ ! -z "$RAWFASTQFILE" -a -f "$RAWFASTQFILE" ]; then
                PER=$(echo $CLIPPED_READS_COUNT | perl -ane 'chomp($_); $per=($_*100)/'$RAW_READS_COUNT'; printf("%0.2f%", $per);')
            fi
        fi

        ## determine unique id
        ID=$(zless $MAPSTATFILE | grep "Map for" | sed 's/Map for //g' | sed 's/\..*//g')

        ## tabulate mapping statistics (header)
        echo -ne "id\t#reads (raw)\t#reads (after qualityCheck)\t#reads (for mapping)\t#reads (unpaired)\t#reads (unmapped)\t#reads (aligned 1 time)\t#reads (aligned >1 time)\talignment rate"
        if [ ! -z "$BAMFILE" ]; then
            ## add header info, if spike-in genome is present
            if [ "$(samtools idxstats $BAMFILE 2>/dev/null | cut -f 1 | grep "_" | wc -l)" -gt 0 ]; then
               samtools idxstats $BAMFILE 2>/dev/null | perl -ane '
                if($F[0]!~/\*/) {
                    $F[0]=~s/^.*\_//g; $count{$F[0]}+=$F[2];
                }
                END {
                    foreach(keys(%count)) {
                        printf("\t%s (mapped)", $_);
                    }
                }'
            fi

            echo -ne "\t#reads (QC-passed)\t#reads (mapped)\t#reads (PCR duplicates)"

            BAMID=$(echo $BAMFILE | sed 's/\.bam$//g')
            if [ "$(ls $ID* | grep "_dupRemoved_" | grep "bam$" | wc -l)" -gt 0 ]; then
                echo $(ls $ID* | grep "_dupRemoved_" | grep "bam$") | perl -ane 'foreach(@F) { $_=~s/^.*\_//g; $_=~s/\..*//g; print "\t$_ (final mapped)"; }'
            elif [ "$(ls $ID* | grep "_dupRemoved" | grep "bam$" | wc -l)" -gt 0 ]; then
                echo $(ls $ID* | grep "_dupRemoved" | grep "bam$") | perl -ane 'foreach(@F) { $_=~s/^.*\_//g; $_=~s/\..*//g; print "\t$_ (final mapped)"; }'
            fi
        fi
        echo
        echo -ne "$ID"

        ## tabulate mapping statistics (values)
        zless $MAPSTATFILE | perl -ane 'BEGIN { print "\t'$RAW_READS_COUNT'\t'$CLIPPED_READS_COUNT' ('$PER')"; } if($_=~/^[0-9\s]+/) { $_=~s/\;.*//g; chomp($_); $_=~s/\s+[a-zA-Z]+.*//g; print "\t$_"; }'
        TOTAL_READS=$(zless $MAPSTATFILE | grep "reads; of these:" | cut -f 1 -d " ")

        if [ ! -z "$BAMFILE" ]; then
            ## add mapping info, if spike-in genome is present
            if [ "$(samtools idxstats $BAMFILE 2>/dev/null | cut -f 1 | grep "_" | wc -l)" -gt 0 ]; then
               samtools idxstats $BAMFILE 2>/dev/null | perl -ane '
                if($F[0]!~/\*/) {
                    $F[0]=~s/^.*\_//g; $count{$F[0]}+=$F[2];
                }
                END {
                    foreach(keys(%count)) {
                        print "\t$count{$_}";
                    }
                }'
            fi
            samtools flagstat $BAMFILE | perl -ane '
                if($_=~/QC-passed/) { $_=~s/\s+.*//g; $qcPassed_reads=$_; }
                elsif($_=~/mapped\s+\(/) { $_=~s/\s+.*//g; $mapped_reads=$_; }
                elsif($_=~/duplicates/) { $_=~s/\s+.*//g; $duplicate_reads=$_; }
                END {
                    printf("\t%s\t%s (%0.2f)\t%s (%0.2f)", $qcPassed_reads, $mapped_reads, ($mapped_reads*100)/$qcPassed_reads, $duplicate_reads, ($duplicate_reads*100)/$qcPassed_reads);
                }'

            BAMID=$(echo $BAMFILE | sed 's/\.bam$//g')
            if [ "$(ls $BAMID* | grep "_dupRemoved_" | grep "bam$" | wc -l)" -gt 0 ]; then
                for i in $(ls $BAMID* | grep "_dupRemoved_" | grep "bam$"); do
                    samtools flagstat $i | perl -ane 'if($_=~/mapped\s+\(/) { $_=~s/\s+.*//g; print "\t$_"; } '
                done
            elif [ "$(ls $BAMID* | grep "_dupRemoved" | grep "bam$" | wc -l)" -gt 0 ]; then
                for i in $(ls $BAMID* | grep "_dupRemoved" | grep "bam$"); do
                    samtools flagstat $i | perl -ane 'if($_=~/mapped\s+\(/) { $_=~s/\s+.*//g; print "\t$_"; } '
                done
            fi

            #samtools view -f 4 $BAMFILE  | cut -f 10 | sort | uniq -c | sort -nr > $ID.unmapped;
        fi
        echo
    ## compute mapping statistics for paired-end data
    else
        ## determine unique id
        ID=$(zless $MAPSTATFILE | grep "Map for" | sed 's/Map for //g' | sed 's/\..*//g')

        ## tabulate mapping statistics (header)
        echo -ne "id\t#reads (for mapping)\t#reads (paired)\t#reads (unmapped)\t#reads (aligned 1 time)\t#reads (aligned >1 time)"
        if [ ! -z "$BAMFILE" ]; then
            ## add header info, if spike-in genome is present
            if [ "$(samtools idxstats $BAMFILE 2>/dev/null | cut -f 1 | grep "_" | wc -l)" -gt 0 ]; then
               samtools idxstats $BAMFILE 2>/dev/null | perl -ane '
                if($F[0]!~/\*/) {
                    $F[0]=~s/^.*\_//g; $count{$F[0]}+=$F[2];
                }
                END {
                    foreach(keys(%count)) {
                        printf("\t%s (mapped)", $_);
                    }
                }'
            fi

            echo -ne "\t#reads (QC-passed)\t#reads (mapped)\t#reads (paired)\t#reads (singleton)\t#reads (PCR duplicates)"

            BAMID=$(echo $BAMFILE | sed 's/\.bam$//g')
            if [ "$(ls $BAMID* | grep "_dupRemoved_" | grep "bam$" | wc -l)" -gt 0 ]; then
                echo $(ls $BAMID* | grep "_dupRemoved_" | grep "bam$") | perl -ane 'foreach(@F) { $_=~s/^.*\_//g; $_=~s/\..*//g; print "\t$_ (final mapped)"; }'
            elif [ "$(ls $BAMID* | grep "_dupRemoved" | grep "bam$" | wc -l)" -gt 0 ]; then
                echo $(ls $BAMID* | grep "_dupRemoved" | grep "bam$") | perl -ane 'foreach(@F) { print "\tgenome (final mapped)"; }'
            fi
        fi
        echo
        echo -ne "$ID"
        
        ## tabulate mapping statistics (values)
        zless $MAPSTATFILE | perl -ane 'if($_=~/\----/) { last; } if($_=~/^[0-9\s]+/) { $_=~s/\;.*//g; chomp($_); $_=~s/\s+[a-zA-Z]+.*//g; print "\t$_"; }'
        TOTAL_READS=$(zless $MAPSTATFILE | grep "reads; of these:" | cut -f 1 -d " ")

        if [ ! -z "$BAMFILE" ]; then
            ## add mapping info, if spike-in genome is present
            if [ "$(samtools idxstats $BAMFILE 2>/dev/null | cut -f 1 | grep "_" | wc -l)" -gt 0 ]; then
               samtools idxstats $BAMFILE 2>/dev/null | perl -ane '
                if($F[0]!~/\*/) {
                    $F[0]=~s/^.*\_//g; $count{$F[0]}+=$F[2];
                }
                END {
                    foreach(keys(%count)) {
                        printf("\t%0.0f", $count{$_}/2);
                    }
                }'
            fi

            samtools flagstat $BAMFILE | perl -ane '
                if($_=~/QC-passed/) { $_=~s/\s+.*//g; $qcPassed_reads=sprintf("%0.0f", $_/2); }
                elsif($_=~/paired in sequencing/) { $_=~s/\s+.*//g; $mapped_reads=sprintf("%0.0f", $_/2); }
                elsif($_=~/properly paired/) { $_=~s/\s+.*//g; $paired_reads=sprintf("%0.0f", $_/2); }
                elsif($_=~/singletons/) { $_=~s/\s+.*//g; $singleton_reads=sprintf("%0.0f", $_/2); }
                elsif($_=~/duplicates/) { $_=~s/\s+.*//g; $duplicate_reads=sprintf("%0.0f", $_/2); }
                END {
                    printf("\t%s\t%s (%0.2f)\t%s (%0.2f)\t%s (%0.2f)\t %s (%0.2f)", $qcPassed_reads, $mapped_reads, ($mapped_reads*100)/$qcPassed_reads, $paired_reads, ($paired_reads*100)/$qcPassed_reads, $singleton_reads, ($singleton_reads*100)/$qcPassed_reads, $duplicate_reads, ($duplicate_reads*100)/$qcPassed_reads);
                }'

            BAMID=$(echo $BAMFILE | sed 's/\.bam$//g')
            if [ "$(ls $BAMID* | grep "_dupRemoved_" | grep "bam$" | wc -l)" -gt 0 ]; then
                for i in $(ls $BAMID* | grep "_dupRemoved_" | grep "bam$"); do
                    samtools flagstat $i | perl -ane 'if($_=~/properly paired/) { $_=~s/\s+.*//g; printf("\t%0.0f", $_/2); } '
                done
            elif [ "$(ls $BAMID* | grep "_dupRemoved" | grep "bam$" | wc -l)" -gt 0 ]; then
                for i in $(ls $BAMID* | grep "_dupRemoved" | grep "bam$"); do
                    samtools flagstat $i | perl -ane 'if($_=~/properly paired/) { $_=~s/\s+.*//g; printf("\t%0.0f", $_/2); } '
                done
            fi

        fi
        echo
    fi
else
    ## determine unique id
    ID=$(zless $MAPSTATFILE | grep "Map for" | sed 's/Map for //g' | sed 's/\..*//g')

    ## tabulate mapping statistics (header)
    echo -e "id\t#reads (raw)\t#reads (aligned 1 time)\t#reads (aligned >1 time)"
    
    ## tabulate mapping statistics (values)
    TOTAL_READS=$(zless $MAPSTATFILE | grep "Number of input reads" | perl -ane 'print $F[scalar(@F)-1]."\n";')

    UNIQUELY_MAPPED=$(zless $MAPSTATFILE | grep "Uniquely mapped reads number" | perl -ane 'print $F[scalar(@F)-1]."\n";')
    UNIQUELY_MAPPED_PER=$(zless $MAPSTATFILE | grep "Uniquely mapped reads %" | perl -ane 'print $F[scalar(@F)-1]."\n";')

    MULTI_MAPPED=$(zless $MAPSTATFILE | grep 'Number of reads mapped to multiple loci\|Number of reads mapped to too many loci' | perl -ane '$sum+=$F[scalar(@F)-1]; END { print "$sum\n"; }')
    MULTI_MAPPED_PER=$(zless $MAPSTATFILE | grep '% of reads mapped to multiple loci\|% of reads mapped to too many loci' | perl -ane '$per=$F[scalar(@F)-1]; $per=~s/\%//g; $sum+=$per; END { print "$sum%\n"; }')

    echo -e "$ID\t$TOTAL_READS\t$UNIQUELY_MAPPED ($UNIQUELY_MAPPED_PER)\t$MULTI_MAPPED ($MULTI_MAPPED_PER)"
fi
exit
