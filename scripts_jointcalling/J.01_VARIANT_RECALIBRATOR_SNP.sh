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
	DBSNP=$5
	HAPMAP=$6
	OMNI_1KG=$7
	HI_CONF_1KG_PHASE1_SNP=$8
	SAMPLE_SHEET=$9
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=${10}
	SEND_TO=${11}

START_VARIANT_RECALIBRATOR_SNP=`date '+%s'` # capture time process starts for wall clock tracking purposes.

	# explicitly state the maximum number of gaussians to start off with

		MAX_GAUSSIANS="8"

	# construct command line

		CMD="singularity exec $GATK_3_7_0_CONTAINER java -jar" \
			CMD=$CMD" /usr/GenomeAnalysisTK.jar" \
		CMD=$CMD" -T VariantRecalibrator" \
			CMD=$CMD" -R $REF_GENOME" \
			CMD=$CMD" --input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf" \
			CMD=$CMD" -recalFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.recal" \
			CMD=$CMD" -tranchesFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.tranches" \
			CMD=$CMD" -rscriptFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.R" \
			CMD=$CMD" -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP" \
			CMD=$CMD" -resource:omni,known=false,training=true,truth=true,prior=12.0 $OMNI_1KG" \
			CMD=$CMD" -resource:1000G,known=false,training=true,truth=false,prior=10.0 $HI_CONF_1KG_PHASE1_SNP" \
			CMD=$CMD" -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP" \
			CMD=$CMD" -mode SNP" \
			CMD=$CMD" --disable_auto_index_creation_and_locking_when_reading_rods" \
			CMD=$CMD" -an QD" \
			CMD=$CMD" -an MQ" \
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
			CMD=$CMD" -tranche 90.0" \
			CMD=$CMD' --maxGaussians '$MAX_GAUSSIANS

		echo $CMD | bash

	# check the exit signal at this point.

		SCRIPT_STATUS=`echo $?`

	# if vqsr fails then retry by decrementing the number of max gaussians by 1 until you get to 4
	# if it still does not work after setting it to four then stop trying

		if [ $SCRIPT_STATUS -ne 0 ]
			then
				until [[ $SCRIPT_STATUS -eq 0 || $MAX_GAUSSIANS -le 1 ]]
					do
						CMD=$(echo $CMD | sed 's/ --maxGaussians '"$MAX_GAUSSIANS"'//g')
						MAX_GAUSSIANS=$[$MAX_GAUSSIANS-1]
						CMD=$CMD' --maxGaussians '$MAX_GAUSSIANS
						echo $CMD | bash
						SCRIPT_STATUS=`echo $?`
				done
		fi

	# if it fails the first time but ultimately works send a notification saying that the parameter has changed and that the methods document needs to change for release
	# if it ends up failing altogether send a notification saying that I need to look at it.
	# will probably have to start with removing the MQ annotation and go from there.

		if [[ $SCRIPT_STATUS -eq 0 && $MAX_GAUSSIANS -ge 1 && $MAX_GAUSSIANS -lt 7 ]]
			then
				printf "The number of max Gaussians has been changed to $MAX_GAUSSIANS for\n \
				PROJECT:\n \
				$PROJECT\n" \
				| mail -s "SNP VariantRecalibrator parameter changed for $PROJECT" \
				$SEND_TO
			elif [ $SCRIPT_STATUS -ne 0 ]
				then
					printf "This has failed SNP VariantRecalibrator and Kurt needs to look at this for:\n \
					PROJECT:\n \
					$PROJECT\n" \
					| mail -s "SNP VariantRecalibrator FAILED for $PROJECT" \
					$SEND_TO
			else
			:
		fi

	# send final working command line to command line file

		echo $CMD >> $CORE_PATH/$PROJECT/COMMAND_LINES/$PROJECT"_command_lines.txt"
		echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$PROJECT"_command_lines.txt"

	# if exit does not equal 0 then exit with whatever the exit signal is at the end.
	# also write to file that this job failed

		if [ "$SCRIPT_STATUS" -ne 0 ]
		 then
			echo $SM_TAG $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
			>> $CORE_PATH/$PROJECT/TEMP/$SAMPLE_SHEET_NAME"_"$SUBMIT_STAMP"_ERRORS.txt"
			exit $SCRIPT_STATUS
		fi

END_VARIANT_RECALIBRATOR_SNP=`date '+%s'` # capture time process ends for wall clock tracking purposes.

# write out timing metrics to file

	echo $PROJECT",J.001,VARIANT_RECALIBRATOR_SNP,"$HOSTNAME","$START_VARIANT_RECALIBRATOR_SNP","$END_VARIANT_RECALIBRATOR_SNP \
	>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

# MOVE TRANCHES PDF FILE TO TEMP TO BE DELETED AFTER A SUCCESSFUL RUN

	mv -v $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.tranches.pdf \
	$CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.HC.SNV.tranches.pdf

# exit with the signal from the program

	exit $SCRIPT_STATUS
