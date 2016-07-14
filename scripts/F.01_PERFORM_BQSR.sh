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

RIS_ID=${SM_TAG%@*}
BARCODE_2D=${SM_TAG#*@}

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
SM_TAG=$5
REF_GENOME=$6
KNOWN_INDEL_1=$7
KNOWN_INDEL_2=$8
DBSNP=$9
BAIT_BED=${10}

## --BQSR using data only from the baited intervals
## --I am actually going to downsample here, b/c it actually makes more sense to do so.

START_PERFORM_BQSR=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $CORE_PATH/$PROJECT/TEMP/$SM_TAG".realign.bam" \
-R $REF_GENOME \
-knownSites $KNOWN_INDEL_1 \
-knownSites $KNOWN_INDEL_2 \
-knownSites $DBSNP \
-L $BAIT_BED \
-nct 8 \
-o $CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT/$SM_TAG"_PERFORM_BQSR.bqsr"

END_PERFORM_BQSR=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",F.001,PERFORM_BQSR,"$HOSTNAME","$START_PERFORM_BQSR","$END_PERFORM_BQSR \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

md5sum $CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT/$SM_TAG"_PERFORM_BQSR.bqsr" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-I $CORE_PATH/$PROJECT/TEMP/$SM_TAG".realign.bam" \
-R $REF_GENOME \
-knownSites $KNOWN_INDEL_1 \
-knownSites $KNOWN_INDEL_2 \
-knownSites $DBSNP \
-L $BAIT_BED \
-nct 8 \
-o $CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT/$SM_TAG"_PERFORM_BQSR.bqsr" \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"
