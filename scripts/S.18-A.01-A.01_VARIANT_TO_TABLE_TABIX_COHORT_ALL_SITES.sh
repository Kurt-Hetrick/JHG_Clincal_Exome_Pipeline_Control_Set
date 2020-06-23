# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q bigdata.q,lemon.q,prod.q,rnd.q,uhoh.q,cgc.q

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

TABIX_DIR=$1
CORE_PATH=$2

PROJECT=$3

# Filter to just on all of the variants all

START_VARIANT_TO_TABLE_TABIX_COHORT=`date '+%s'`

# not doing --splitMultiallelic here...maybe do one as an example and discuss with Molly
# do an example of molten output to look at/show molly

$TABIX_DIR/tabix -s 1 -b 2 -e 2 -c C \
$CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_SET.VQSR.ANNOTATED.txt.gz

END_VARIANT_TO_TABLE_TABIX_COHORT=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",V.01,VARIANT_TO_TABLE_CONTROL_SET_ALL_SITES,"$HOSTNAME","$START_VARIANT_TO_TABLE_TABIX_COHORT","$START_VARIANT_TO_TABLE_TABIX_COHORT \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $TABIX_DIR/tabix -s 1 -b 2 -e 2 -c C \
$CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_SET.VQSR.ANNOTATED.txt.gz \
>> $CORE_PATH/$PROJECT/$PROJECT".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$PROJECT".COMMAND.LINES.txt"

md5sum $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_SET.VQSR.ANNOTATED.txt.gz.tbi \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
