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

	GATK_3_7_0_CONTAINER=$1
	CORE_PATH=$2

	PROJECT=$3
	REF_GENOME=$4
	MILLS_1KG_GOLD_INDEL=$5
	$SAMPLE_SHEET=$6
	$SUBMIT_STAMP=$7

START_VARIANT_RECALIBRATOR_INDEL=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# construct command line

		CMD="singularity exec $GATK_3_7_0_CONTAINER java -jar" \
			CMD=$CMD" /usr/GenomeAnalysisTK.jar" \
		CMD=$CMD" -T VariantRecalibrator" \
			CMD=$CMD" -R $REF_GENOME" \
			CMD=$CMD" --input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf" \
			CMD=$CMD" -recalFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.INDEL.recal" \
			CMD=$CMD" -tranchesFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.INDEL.tranches" \
			CMD=$CMD" -rscriptFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.INDEL.R" \
			CMD=$CMD" -resource:mills,known=true,training=true,truth=true,prior=12.0 $MILLS_1KG_GOLD_INDEL" \
			CMD=$CMD" -mode INDEL" \
			CMD=$CMD" --maxGaussians 4" \
			CMD=$CMD" --disable_auto_index_creation_and_locking_when_reading_rods" \
			CMD=$CMD" -an QD" \
			CMD=$CMD" -an MQRankSum" \
			CMD=$CMD" -an ReadPosRankSum" \
			CMD=$CMD" -an FS" \
			CMD=$CMD" -an SOR" \
			CMD=$CMD" -tranche 100.0" \
			CMD=$CMD" -tranche 99.9" \
			CMD=$CMD" -tranche 99.8" \
			CMD=$CMD" -tranche 99.7" \
			CMD=$CMD" -tranche 99.6" \
			CMD=$CMD" -tranche 99.5" \
			CMD=$CMD" -tranche 99.4" \
			CMD=$CMD" -tranche 99.3" \
			CMD=$CMD" -tranche 99.2" \
			CMD=$CMD" -tranche 99.1" \
			CMD=$CMD" -tranche 99.0" \
			CMD=$CMD" -tranche 98.0" \
			CMD=$CMD" -tranche 97.0" \
			CMD=$CMD" -tranche 96.0" \
			CMD=$CMD" -tranche 95.0" \
			CMD=$CMD" -tranche 90.0"

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

END_VARIANT_RECALIBRATOR_INDEL=`date '+%s'` # capture time process starts for wall clock tracking purposes.

# write out timing metrics to file

	echo $PROJECT",J.001,VARIANT_RECALIBRATOR_INDEL,"$HOSTNAME","$START_VARIANT_RECALIBRATOR_INDEL","$END_VARIANT_RECALIBRATOR_INDEL \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# exit with the signal from the program

	exit $SCRIPT_STATUS
