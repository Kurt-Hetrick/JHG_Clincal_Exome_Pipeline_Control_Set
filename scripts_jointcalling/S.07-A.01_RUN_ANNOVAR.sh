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

	CIDRSEQSUITE_ANNOVAR_JAVA=$1
	CIDRSEQSUITE_DIR=$2
	CIDRSEQSUITE_PROPS_DIR=$3
	CORE_PATH=$4

	PROJECT=$5
	SM_TAG=$6
	SAMPLE_SHEET=$7
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$8

# Copy the vcf that needs to have annovar run on it to the sample specific temp folder.

	cp -rvf $CORE_PATH/$PROJECT/VCF/FILTERED_ON_BAIT/$SM_TAG".VARIANT_SITES.vcf" \
	$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"

# Run annovar like a boss

START_RUN_ANNOVAR=`date '+%s'` # capture time process starts for wall clock tracking purposes

	# construct command line

		CMD="$CIDRSEQSUITE_ANNOVAR_JAVA/java -jar -Duser.home=$CIDRSEQSUITE_PROPS_DIR" \
		CMD=$CMD" $CIDRSEQSUITE_DIR/CIDRSeqSuite.jar" \
			CMD=$CMD" -pipeline" \
			CMD=$CMD" -annovar_directory_annotation" \
			CMD=$CMD" $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"" \
			CMD=$CMD" $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR""

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

END_RUN_ANNOVAR=`date '+%s'` # capture time process ends for wall clock tracking

# write out timing metrics to file

	echo $PROJECT",Q.001,RUN_ANNOVAR,"$HOSTNAME","$START_RUN_ANNOVAR","$END_RUN_ANNOVAR \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
