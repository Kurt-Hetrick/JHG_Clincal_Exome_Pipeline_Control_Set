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

	ALIGNMENT_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	REF_GENOME=$4
	SAMPLE_SHEET=$5
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$6

# Filter to just on passing variants for all samples

START_FILTER_ALL_CONTROLS_VARIANT_ONLY_PASS=`date '+%s'` # capture time process starts for wall clock tracking purposes

	# construct command line

		CMD="singularity exec $ALIGNMENT_CONTAINER java -jar" \
			CMD=$CMD" /gatk/gatk.jar" \
		CMD=$CMD" SelectVariants" \
			CMD=$CMD" --reference $REF_GENOME" \
			CMD=$CMD" --variant $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz" \
			CMD=$CMD" --output $CORE_PATH/$PROJECT/JOINT_VCF/CONTROLS_DATA_SET.VQSR.ANNOTATED.VARIANT_ONLY_PASS.vcf.gz" \
			CMD=$CMD" --exclude-non-variants" \
			CMD=$CMD" --exclude-filtered"

	# write command line to file and execute the command line

		echo $CMD >> $CORE_PATH/$PROJECT/COMMAND_LINES/$PROJECT"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$PROJECT"_command_lines.txt"
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

END_FILTER_ALL_CONTROLS_VARIANT_ONLY_PASS=`date '+%s'` # capture time process ends for wall clock tracking purposes.

# write out timing metrics to file

	echo $PROJECT",S.001,FILTER_COHORT_VARIANT_ONLY_PASS,"$HOSTNAME","$START_FILTER_COHORT_VARIANT_ONLY_PASS","$END_FILTER_COHORT_VARIANT_ONLY_PASS \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
