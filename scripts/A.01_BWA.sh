# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash 

# tell sge to submit any of these queue when available
#$ -q rnd.q,prod.q,test.q

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

BWA_DIR=$1
JAVA_1_8=$2
PICARD_DIR=$3
CORE_PATH=$4

PROJECT=$5 # the Seq Proj folder name. 1st column in sample manifest
FLOWCELL=$6 # flowcell that sample read group was performed on. 2nd column of sample manifest
LANE=$7 # lane of flowcell that sample read group was performed on. 3rd column of the sample manifest
INDEX=$8 # sample barcode. 4th column of the sample manifest
PLATFORM=$9 # type of sequencing chemistry matching SAM specification. 5th column of the sample manifest.
LIBRARY_NAME=${10} # library group of the sample read group.
								# Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not
								# 6th column of the sample manifest
RUN_DATE=${11} # should be the run set up date to match the seq run folder name, but it has been arbitrarily populated. field X of manifest.
SM_TAG=${12} # sample ID. sample name for all files, etc. field X of manifest
CENTER=${13} # the center/funding mechanism. field X of manifest.
DESCRIPTION=${14} # Generally we use to denote the sequencer setting (e.g. rapid run). field X of manifest.
REF_GENOME=${15} # the reference genome used in the analysis pipeline. field X of manifest.

RIS_ID=${SM_TAG%@*} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
BARCODE_2D=${SM_TAG#*@} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.

PLATFORM_UNIT=$FLOWCELL"_"$LANE"_"$INDEX

# Need to convert data in sample manifest to Iso 8601 date since we are not using bwa mem to populate this.
# Picard AddOrReplaceReadGroups is much more stringent here.

ISO_8601=`echo $RUN_DATE \
|awk '{split ($0,DATES,"/"); \
if (length(DATES[1]) < 2 && length(DATES[2]) < 2) \
print DATES[3]"-0"DATES[1]"-0"DATES[2]"T00:00:00-0500"; \
else if (length(DATES[1]) < 2 && length(DATES[2]) > 1) \
print DATES[3]"-0"DATES[1]"-"DATES[2]"T00:00:00-0500"; \
else if(length(DATES[1]) > 1 && length(DATES[2]) < 2) \
print DATES[3]"-"DATES[1]"-0"DATES[2]"T00:00:00-0500"; \
else print DATES[3]"-"DATES[1]"-"DATES[2]"T00:00:00-0500"}'`

# -----Alignment and BAM post-processing-----

# --bwa mem
# --pipe to MergeSamFiles to sort and write a bam file.--

# look for fastq files. allow fastq.gz and fastq extensions.

FASTQ_1=`ls $CORE_PATH/$PROJECT/FASTQ/$PLATFORM_UNIT"_1.fastq"*`
FASTQ_2=`ls $CORE_PATH/$PROJECT/FASTQ/$PLATFORM_UNIT"_2.fastq"*`

# BWA POPULATES SEQUENCE DICTIONARY...MIGHT CONSIDER FILLING THIS MORE COMPLETELY...LOW PRIORITY

START_BWA_MEM=`date '+%s'`

$BWA_DIR/bwa mem \
-M \
-t 4 \
$REF_GENOME \
$FASTQ_1 \
$FASTQ_2 \
| $JAVA_1_8/java -jar \
$PICARD_DIR/picard.jar AddOrReplaceReadGroups \
INPUT=/dev/stdin \
CREATE_INDEX=true \
SORT_ORDER=coordinate \
RGID= $FLOWCELL"_"$LANE \
RGLB= $LIBRARY_NAME \
RGPL= $PLATFORM \
RGPU= $PLATFORM_UNIT \
RGSM= $SM_TAG \
RGCN= $CENTER \
RGDS= $DESCRIPTION \
RGDT= $ISO_8601 \
RGPG= CGC_CIDRSEQSUITE_0.0.0-0 \
OUTPUT=$CORE_PATH/$PROJECT/TEMP/$PLATFORM_UNIT".bam"

END_BWA_MEM=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",A.001,BWA_MEM,"$HOSTNAME","$START_BWA_MEM","$END_BWA_MEM \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# I'm guessing the pipe screws up the echo. Need to look up how to echo a pipe.

echo $BWA_DIR/bwa mem \
-M \
-t 4 \
$REF_GENOME \
$FASTQ_1 \
$FASTQ_2 \
\| $JAVA_1_8/java -jar \
$PICARD_DIR/picard.jar AddOrReplaceReadGroups \
INPUT=/dev/stdin \
CREATE_INDEX=true \
SORT_ORDER=coordinate \
RGID= $FLOWCELL"_"$LANE \
RGLB= $LIBRARY_NAME \
RGPL= $PLATFORM \
RGPU= $PLATFORM_UNIT \
RGSM= $SM_TAG \
RGCN= $CENTER \
RGDS= $DESCRIPTION \
RGDT= $ISO_8601 \
RGPG= CGC_CIDRSEQSUITE_0.0.0-0 \
OUTPUT=$CORE_PATH/$PROJECT/TEMP/$PLATFORM_UNIT".bam" \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# RGPG # PROGRAM...HMM...CIDRSEQSUITE VERSION MAYBE?...# Actually this is a really good idea for a clinical pipeline
# RGPM # THIS IS WHAT WE ACTUALLY USE DESCRIPTION FOR...
