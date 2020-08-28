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

	BWA_DIR=$1
	SAMBLASTER_DIR=$2
	JAVA_1_8=$3
	PICARD_DIR=$4

	CORE_PATH=$5
	PROJECT=$6
	FLOWCELL=$7
	LANE=$8
	INDEX=$9
		PLATFORM_UNIT=$FLOWCELL"_"$LANE"_"$INDEX
		FIXED_PLATFORM_UNIT=`echo $PLATFORM_UNIT | sed 's/~/*/g'`
	PLATFORM=${10}
	LIBRARY_NAME=${11}
	RUN_DATE=${12}
	SM_TAG=${13}
	CENTER=${14}
	SEQUENCER_MODEL=${15}
	REF_GENOME=${16}
	PIPELINE_VERSION=${17}
	BAIT_BED=${18}
		BAIT_NAME=$(basename $BAIT_BED .bed)
	TARGET_BED=${19}
		TARGET_NAME=$(basename $TARGET_BED .bed)
	TITV_BED=${20}
		TITV_NAME=$(basename $TITV_BED .bed)
	SAMPLE_SHEET=${21}
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=${22}
	NOVASEQ_REPO=${23}

# Need to convert data in sample manifest to Iso 8601 date since we are not using bwa mem to populate this.
# Picard AddOrReplaceReadGroups is much more stringent here.

	if [[ $RUN_DATE = *"-"* ]];
		then

			# for when the date is this 2018-09-05

				ISO_8601=`echo $RUN_DATE \
					| awk '{print "'$RUN_DATE'" "T00:00:00-0500"}'`

		else

			# for when the data is like this 4/26/2018

				ISO_8601=`echo $RUN_DATE \
					| awk '{split ($0,DATES,"/"); \
					if (length(DATES[1]) < 2 && length(DATES[2]) < 2) \
					print DATES[3]"-0"DATES[1]"-0"DATES[2]"T00:00:00-0500"; \
					else if (length(DATES[1]) < 2 && length(DATES[2]) > 1) \
					print DATES[3]"-0"DATES[1]"-"DATES[2]"T00:00:00-0500"; \
					else if(length(DATES[1]) > 1 && length(DATES[2]) < 2) \
					print DATES[3]"-"DATES[1]"-0"DATES[2]"T00:00:00-0500"; \
					else print DATES[3]"-"DATES[1]"-"DATES[2]"T00:00:00-0500"}'`
	fi

# look for fastq files. allow fastq.gz and fastq extensions.
# If NovaSeq is contained in the Description field in the sample sheet then assume that ILMN BCL2FASTQ is used.
# Files are supposed to be in /mnt/instrument_files/novaseq/Run_Folder/FASTQ/Project/
# FILENAME-> 137233-0238091146_S49_L002_R1_001.fastq.gz	(SMTAG_ASampleIndexOfSomeSort_4DigitLane_Read_literally001.fastq.gz)
# Otherwise assume that files are demultiplexed with cidrseqsuite and follow previous naming conventions.
# I got files from yale, that used the illumina naming conventions and actually went a step farther and broke files by tile (i think).
	## I concatenated them and then added 000 for the tile so added that to the end of the non novaseq fastq file look up

	if [[ $SEQUENCER_MODEL == *"NovaSeq"* ]]
		then

			NOVASEQ_RUN_FOLDER=`ls $NOVASEQ_REPO | grep $FLOWCELL`

			FINDPATH=$NOVASEQ_REPO/$NOVASEQ_RUN_FOLDER/FASTQ/$PROJECT

			# look for illumina file naming convention for novaseq flowcells
			# if it is found in the project/fastq folder under active, then use that one
			FASTQ_1=`( echo du --max-depth=1 -a $FINDPATH/$SM_TAG"*" -a $FINDPATH/$FIXED_PLATFORM_UNIT"*" 2\> /dev/null \| grep "L00"$LANE"_R1_001.fastq" \| cut -f 2 | bash ; \
				ls $CORE_PATH/$PROJECT/FASTQ/$FIXED_PLATFORM_UNIT"_1.fastq"* 2> /dev/null) | tail -n 1`
			FASTQ_2=`( echo du --max-depth=1 -a $FINDPATH/$SM_TAG"*" -a $FINDPATH/$FIXED_PLATFORM_UNIT"*" 2\> /dev/null \| grep "L00"$LANE"_R2_001.fastq" \| cut -f 2 | bash ; \
				ls $CORE_PATH/$PROJECT/FASTQ/$FIXED_PLATFORM_UNIT"_2.fastq"* 2> /dev/null) | tail -n 1`

		else
			FASTQ_1=`(ls $CORE_PATH/$PROJECT/FASTQ/$FIXED_PLATFORM_UNIT"_1.fastq"* 2> /dev/null ; ls $CORE_PATH/$PROJECT/FASTQ/$FIXED_PLATFORM_UNIT"_R1_000.fastq"* 2> /dev/null; ls $CORE_PATH/$PROJECT/FASTQ/$SM_TAG"_R1_001.fastq"* 2> /dev/null; ls $CORE_PATH/$PROJECT/FASTQ/$SM_TAG"_1.fastq"* 2> /dev/null)`
			FASTQ_2=`(ls $CORE_PATH/$PROJECT/FASTQ/$FIXED_PLATFORM_UNIT"_2.fastq"* 2> /dev/null ; ls $CORE_PATH/$PROJECT/FASTQ/$FIXED_PLATFORM_UNIT"_R2_000.fastq"* 2> /dev/null; ls $CORE_PATH/$PROJECT/FASTQ/$SM_TAG"_R2_001.fastq"* 2> /dev/null; ls $CORE_PATH/$PROJECT/FASTQ/$SM_TAG"_2.fastq"* 2> /dev/null)`
	fi

# -----Alignment and BAM post-processing-----

	# bwa mem
	# pipe to samblaster to add MC, etc tags
	# pipe to AddOrReplaceReadGroups to populate the header--

# bwa mem for paired end reads

	BWA_PE ()
		{
			START_BWA_MEM=`date '+%s'`
				# if any part of pipe fails set exit to non-zero

				set -o pipefail

				$BWA_DIR/bwa mem \
					-K 100000000 \
					-Y \
					-t 4 \
					$REF_GENOME \
					$FASTQ_1 \
					$FASTQ_2 \
				| $SAMBLASTER_DIR/samblaster \
					--addMateTags \
					-a \
				| $JAVA_1_8/java -jar \
				$PICARD_DIR/picard.jar \
				AddOrReplaceReadGroups \
				INPUT=/dev/stdin \
				CREATE_INDEX=true \
				SORT_ORDER=queryname \
				RGID=$FLOWCELL"_"$LANE \
				RGLB=$LIBRARY_NAME \
				RGPL=$PLATFORM \
				RGPU=$PLATFORM_UNIT \
				RGPM=$SEQUENCER_MODEL \
				RGSM=$SM_TAG \
				RGCN=$CENTER \
				RGDT=$ISO_8601 \
				RGPG="CIDR_WES-"$PIPELINE_VERSION \
				RGDS=$BAIT_NAME","$TARGET_NAME","$TITV_NAME \
				OUTPUT=$CORE_PATH/$PROJECT/TEMP/$PLATFORM_UNIT".bam"

				# check the exit signal at this point.

					SCRIPT_STATUS=`echo $?`

				# if exit does not equal 0 then exit with whatever the exit signal is at the end.
				# also write to file that this job failed
				# so if it crashes, I just straight out exit
					### ...at first I didn't remember why would I chose that, but I am cool with it
					### ...not good for debugging, but I don't want cmd lines and times when jobs crash tbh if the plan is to possibly distribute them

					if [ "$SCRIPT_STATUS" -ne 0 ]
					 then
						echo $SM_TAG $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
						>> $CORE_PATH/$PROJECT/TEMP/$SAMPLE_SHEET_NAME"_"$SUBMIT_STAMP"_ERRORS.txt"
						exit $SCRIPT_STATUS
					fi

			END_BWA_MEM=`date '+%s'`

			echo $SM_TAG"_"$PROJECT",A.01,BWA_MEM,"$HOSTNAME","$START_BWA_MEM","$END_BWA_MEM \
			>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

			echo $BWA_DIR/bwa mem \
			-K 100000000 \
			-Y \
			-t 4 \
			$REF_GENOME \
			$FASTQ_1 \
			$FASTQ_2 \
			\| $SAMBLASTER_DIR/samblaster \
			--addMateTags \
			-a \
			\| $JAVA_1_8/java -jar \
			$PICARD_DIR/picard.jar \
			AddOrReplaceReadGroups \
			INPUT=/dev/stdin \
			CREATE_INDEX=true \
			SORT_ORDER=queryname \
			RGID=$FLOWCELL"_"$LANE \
			RGLB=$LIBRARY_NAME \
			RGPL=$PLATFORM \
			RGPU=$PLATFORM_UNIT \
			RGPM=$SEQUENCER_MODEL \
			RGSM=$SM_TAG \
			RGCN=$CENTER \
			RGDT=$ISO_8601 \
			RGPG="CIDR_WES-"$PIPELINE_VERSION \
			RGDS=$BAIT_NAME","$TARGET_NAME","$TITV_NAME \
			OUTPUT=$CORE_PATH/$PROJECT/TEMP/$PLATFORM_UNIT".bam" \
			>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

			echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

		}

# bwa mem for single end reads

	BWA_SE ()
		{
			START_BWA_MEM=`date '+%s'`
				# if any part of pipe fails set exit to non-zero

				set -o pipefail

				$BWA_DIR/bwa mem \
					-K 100000000 \
					-Y \
					-t 4 \
					$REF_GENOME \
					$FASTQ_1 \
				| $SAMBLASTER_DIR/samblaster \
					--addMateTags \
					-a \
				| $JAVA_1_8/java -jar \
				$PICARD_DIR/picard.jar \
				AddOrReplaceReadGroups \
				INPUT=/dev/stdin \
				CREATE_INDEX=true \
				SORT_ORDER=queryname \
				RGID=$FLOWCELL"_"$LANE \
				RGLB=$LIBRARY_NAME \
				RGPL=$PLATFORM \
				RGPU=$PLATFORM_UNIT \
				RGPM=$SEQUENCER_MODEL \
				RGSM=$SM_TAG \
				RGCN=$CENTER \
				RGDT=$ISO_8601 \
				RGPG="CIDR_WES-"$PIPELINE_VERSION \
				RGDS=$BAIT_NAME","$TARGET_NAME","$TITV_NAME \
				OUTPUT=$CORE_PATH/$PROJECT/TEMP/$PLATFORM_UNIT".bam"

				# check the exit signal at this point.

					SCRIPT_STATUS=`echo $?`

				# if exit does not equal 0 then exit with whatever the exit signal is at the end.
				# also write to file that this job failed
				# so if it crashes, I just straight out exit
					### ...at first I didn't remember why would I chose that, but I am cool with it
					### ...not good for debugging, but I don't want cmd lines and times when jobs crash tbh if the plan is to possibly distribute them

					if [ "$SCRIPT_STATUS" -ne 0 ]
					 then
						echo $SM_TAG $HOSTNAME $JOB_NAME $USER $SCRIPT_STATUS $SGE_STDERR_PATH \
						>> $CORE_PATH/$PROJECT/TEMP/$SAMPLE_SHEET_NAME"_"$SUBMIT_STAMP"_ERRORS.txt"
						exit $SCRIPT_STATUS
					fi

			END_BWA_MEM=`date '+%s'`

			echo $SM_TAG"_"$PROJECT",A.01,BWA_MEM,"$HOSTNAME","$START_BWA_MEM","$END_BWA_MEM \
			>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

			echo $BWA_DIR/bwa mem \
			-K 100000000 \
			-Y \
			-t 4 \
			$REF_GENOME \
			$FASTQ_1 \
			\| $SAMBLASTER_DIR/samblaster \
			--addMateTags \
			-a \
			\| $JAVA_1_8/java -jar \
			$PICARD_DIR/picard.jar \
			AddOrReplaceReadGroups \
			INPUT=/dev/stdin \
			CREATE_INDEX=true \
			SORT_ORDER=queryname \
			RGID=$FLOWCELL"_"$LANE \
			RGLB=$LIBRARY_NAME \
			RGPL=$PLATFORM \
			RGPU=$PLATFORM_UNIT \
			RGPM=$SEQUENCER_MODEL \
			RGSM=$SM_TAG \
			RGCN=$CENTER \
			RGDT=$ISO_8601 \
			RGPG="CIDR_WES-"$PIPELINE_VERSION \
			RGDS=$BAIT_NAME","$TARGET_NAME","$TITV_NAME \
			OUTPUT=$CORE_PATH/$PROJECT/TEMP/$PLATFORM_UNIT".bam" \
			>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

			echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

		}

# If there is no read 2 the run bwa mem for single end reads, else do paired end.

	if [ -z "$FASTQ_2" ]
		then
		      BWA_SE
		else
		      BWA_PE
	fi


# exit with the signal from the program

	exit $SCRIPT_STATUS
