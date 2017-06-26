# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q cgc.q

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on
# redirecting stderr/stdout to file as a log.

set

CORE_PATH=$1
BEDTOOLS_DIR=$2
CODING_BED=$3

PROJECT=$4
SM_TAG=$5

# PAD THE ANNOTATED CODING BED FILE
# The input bed file could be a variable name based on the padding length

awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10,$4,$5,$6,$7}' $CODING_BED \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOTATED_PADDED_CODING.bed"

START_PER_BASE=`date '+%s'`

awk 'BEGIN {FS=","; OFS="\t"} NR>1 {split($1,BASE,":"); print BASE[1],BASE[2]-1,BASE[2],$2}' \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.EveryBase.csv" \
| $BEDTOOLS_DIR/bedtools intersect \
-a $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOTATED_PADDED_CODING.bed" \
-b - \
-wb \
| awk '{OFS="\t"} BEGIN {print "#CHROM" "\t" "POS" "\t" "GENE" "\t" "TRANSCRIPT" "\t" "EXON" "\t" "STRAND" "\t" "DEPTH"} {print $1,$3,$4,$5,$6,$7,$11}' \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt"

(head -n 1 $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt" ; \
awk '$1~/^[0-9]/' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt" | sort -k1,1n -k 2,2n ; \
awk '$1=="X"' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt" | sort -k1,1n -k 2,2n ; \
awk '$1=="Y"' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt" | sort -k1,1n -k 2,2n ; \
awk '$1=="MT"' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt" | sort -k1,1n -k 2,2n) \
>| $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt"

END_PER_BASE=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.001,REFSEQ_PER_BASE,"$HOSTNAME","$START_PER_CODING","$END_PER_BASE \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# echo awk 'BEGIN {FS=","; OFS="\t"} NR>1 {split($1,BASE,":"); print BASE[1],BASE[2]-1,BASE[2],$2}' \
# $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.csv" \
# | $BEDTOOLS_DIR/bedtools intersect \
# -a $CORE_PATH/$PROJECT/TEMP/PADDED_$CODING_BED \
# -b /dev/stdin \
# -wb \
# | awk '{OFS="\t"} BEGIN {print "#CHROM" "\t" "POS" "\t" "GENE" "\t" "TRANSCIPT" "\t" "EXON" "\t" "STRAND" "\t" "DEPTH"} \
# {print $1,$3,$4,$5,$6,$7,$11}' \
# \>\| $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG"_REFSEQ_PLUS_10bp_PAD.PER.BASE.REPORT.txt" \
# >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
