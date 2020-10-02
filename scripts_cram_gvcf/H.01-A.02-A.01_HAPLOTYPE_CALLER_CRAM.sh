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

	set

	echo

# INPUT VARIABLES

	SAMTOOLS_DIR=$1
	CORE_PATH=$2

	PROJECT=$3
	SM_TAG=$4
	REF_GENOME=$5
	SAMPLE_SHEET=$6
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$7

## --write lossless cram file.

START_HC_CRAM=`date '+%s'`

	$SAMTOOLS_DIR/samtools \
	view \
	-C $CORE_PATH/$PROJECT/TEMP/$SM_TAG".HC.bam" \
	-T $REF_GENOME \
	-@ 4 \
	-o $CORE_PATH/$PROJECT/HC_CRAM/$SM_TAG".HC.cram"

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

	# if exit does not equal 0 then exit with whatever the exit signal is at the end.
	# also write to file that this job failed

		if [ "$SCRIPT_STATUS" -ne 0 ]
		 then
			echo $SM_TAG $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
			>> $CORE_PATH/$PROJECT/TEMP/$SAMPLE_SHEET_NAME"_"$SUBMIT_STAMP"_ERRORS.txt"
			exit $SCRIPT_STATUS
		fi

END_HC_CRAM=`date '+%s'`

echo $SM_TAG"_"$PROJECT",H.01-A.01-A.01,HC_CRAM,"$HOSTNAME","$START_HC_CRAM","$END_HC_CRAM \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $SAMTOOLS_DIR/samtools \
view \
-C $CORE_PATH/$PROJECT/TEMP/$SM_TAG".HC.bam" \
-T $REF_GENOME \
-@ 4 \
-o $CORE_PATH/$PROJECT/HC_CRAM/$SM_TAG".HC.cram" \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# exit with the signal from the program

	exit $SCRIPT_STATUS
