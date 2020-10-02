# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

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

	echo

# INPUT VARIABLES

	SAMTOOLS_DIR=$1
	CORE_PATH=$2

	PROJECT=$3
	SM_TAG=$4
	REF_GENOME=$5

## --index the cram file

START_INDEX_HC_CRAM=`date '+%s'`

	$SAMTOOLS_DIR/samtools \
	index \
	$CORE_PATH/$PROJECT/HC_CRAM/$SM_TAG".HC.cram"

END_INDEX_HC_CRAM=`date '+%s'`

echo $SM_TAG"_"$PROJECT",H.01-A.01-A.01-A.01,INDEX_HC_CRAM,"$HOSTNAME","$START_INDEX_HC_CRAM","$END_INDEX_HC_CRAM \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $SAMTOOLS_DIR/samtools \
index \
$CORE_PATH/$PROJECT/HC_CRAM/$SM_TAG".HC.cram" \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# make a copy/rename the cram index file since their appears to be two useable standards

cp -rvf $CORE_PATH/$PROJECT/HC_CRAM/$SM_TAG".HC.cram.crai" \
$CORE_PATH/$PROJECT/HC_CRAM/$SM_TAG".HC.crai"
