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

	JAVA_1_8=$1
	GATK_DIR_4011=$2
	CORE_PATH=$3

	PROJECT=$4
	SM_TAG=$5
	REF_GENOME=$6
	SAMPLE_SHEET=$7
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$8

## --write out bam file with a 4 bin qscore scheme, remove indel Q scores, emit original Q scores
# have to change the way to specify this jar file eventually. gatk 4 devs are monsters.

START_FINAL_BAM=`date '+%s'`

	$JAVA_1_8/java -jar \
	$GATK_DIR_4011/gatk-package-4.0.1.1-local.jar \
	ApplyBQSR \
	--add-output-sam-program-record \
	--use-original-qualities \
	--emit-original-quals \
	--reference $REF_GENOME \
	--input $CORE_PATH/$PROJECT/TEMP/$SM_TAG".dup.bam" \
	--bqsr-recal-file $CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT/$SM_TAG"_PERFORM_BQSR.bqsr" \
	--static-quantized-quals 10 \
	--static-quantized-quals 20 \
	--static-quantized-quals 30 \
	--output $CORE_PATH/$PROJECT/TEMP/$SM_TAG".bam"

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

END_FINAL_BAM=`date '+%s'`

echo $SM_TAG"_"$PROJECT",E.01,FINAL_BAM,"$HOSTNAME","$START_FINAL_BAM","$END_FINAL_BAM \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar \
$GATK_DIR_4011/gatk-package-4.0.1.1-local.jar \
ApplyBQSR \
--add-output-sam-program-record \
--use-original-qualities \
--emit-original-quals \
--reference $REF_GENOME \
--input $CORE_PATH/$PROJECT/TEMP/$SM_TAG".dup.bam" \
--bqsr-recal-file $CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT/$SM_TAG"_PERFORM_BQSR.bqsr" \
--static-quantized-quals 10 \
--static-quantized-quals 20 \
--static-quantized-quals 30 \
--output $CORE_PATH/$PROJECT/TEMP/$SM_TAG".bam" \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# exit with the signal from the program

	exit $SCRIPT_STATUS
