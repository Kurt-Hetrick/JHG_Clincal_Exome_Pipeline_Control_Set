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

CORE_PATH=$1
VERIFY_DIR=$2

PROJECT=$3
SM_TAG=$4

RIS_ID=${SM_TAG%@*}
BARCODE_2D=${SM_TAG#*@}

## --Running verifyBamID--

START_VERIFYBAMID=`date '+%s'`

$VERIFY_DIR/verifyBamID \
--bam $CORE_PATH/$PROJECT/BAM/$SM_TAG".bam" \
--vcf $CORE_PATH/$PROJECT/TEMP/$SM_TAG".VerifyBamID.vcf" \
--out $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG \
--precise \
--verbose \
--maxDepth 2500

END_VERIFYBAMID=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.001,VERIFYBAMID,"$HOSTNAME","$START_VERIFYBAMID","$END_VERIFYBAMID \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $VERIFY_DIR/verifyBamID \
--bam $CORE_PATH/$PROJECT/BAM/$SM_TAG".bam" \
--vcf $CORE_PATH/$PROJECT/TEMP/$SM_TAG".VerifyBamID.vcf" \
--out $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG \
--precise \
--verbose \
--maxDepth 2500 \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG".selfSM" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG".selfRG" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG".depthSM" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG".depthRG" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID/$SM_TAG".log" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
