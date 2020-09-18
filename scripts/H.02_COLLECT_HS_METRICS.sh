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

	# I'M KEEPING TARGET BED AS AN INPUT VARIABLE EVEN THOUGH IT IS NO LONGER USED EXPLICITLY
	# MORE OF A CUE IN THE WRAPPER SCRIPT AS TO WHAT IS GOING ON

	ALIGNMENT_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	SM_TAG=$4
	REF_GENOME=$5
	BAIT_BED=$6
		BAIT_BED_NAME=(`basename $BAIT_BED .bed`)
	TITV_BED=$7
		TITV_BED_NAME=(`basename $TITV_BED .bed`)
	SAMPLE_SHEET=$8
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$9

# Run Collect HS metrics which generates hybridization metrics for the qc report
## Also generates a per target interval coverage summary

START_COLLECT_HS_METRICS=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec $ALIGNMENT_CONTAINER java -jar" \
			CMD=$CMD" /gatk/picard.jar" \
		CMD=$CMD" CollectHsMetrics" \
		CMD=$CMD" INPUT=$CORE_PATH/$PROJECT/CRAM/$SM_TAG".cram"" \
		CMD=$CMD" OUTPUT=$CORE_PATH/$PROJECT/REPORTS/HYB_SELECTION/$SM_TAG"_hybridization_selection_metrics.txt"" \
		CMD=$CMD" PER_TARGET_COVERAGE=$CORE_PATH/$PROJECT/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE/$SM_TAG"_per_target_coverage.txt"" \
		CMD=$CMD" REFERENCE_SEQUENCE=$REF_GENOME" \
		CMD=$CMD" BAIT_INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME"-picard.bed"" \
		CMD=$CMD" TARGET_INTERVALS=$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$TITV_BED_NAME"-picard.bed"" \
		CMD=$CMD" MINIMUM_MAPPING_QUALITY=20" \
		CMD=$CMD" MINIMUM_BASE_QUALITY=10" \
		CMD=$CMD" BAIT_SET_NAME=$BAIT_BED_NAME" \
		CMD=$CMD" VALIDATION_STRINGENCY=SILENT" \
		CMD=$CMD" COVERAGE_CAP=75000"

	# write command line to file and execute the command line

		echo $CMD >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG"_command_lines.txt"
		echo $CMD | bash

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

END_COLLECT_HS_METRICS=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# write out timing metrics to file

	echo $SM_TAG"_"$PROJECT"_BAM_REPORTS,Z.01,COLLECT_HS_METRICS,"$HOSTNAME","$START_COLLECT_HS_METRICS","$END_COLLECT_HS_METRICS \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
