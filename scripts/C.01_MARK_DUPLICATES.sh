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
	PICARD_DIR=$2
	SAMBAMBA_DIR=$3
	CORE_PATH=$4

	PROJECT=$5
	SM_TAG=$6
	SAMPLE_SHEET=$7
		SAMPLE_SHEET_NAME=$(basename $SAMPLE_SHEET .csv)
	SUBMIT_STAMP=$8
	SEQUENCER_MODEL=$9

	INPUT_BAM_FILE_STRING=${10}
		INPUT=`echo $INPUT_BAM_FILE_STRING | sed 's/,/ /g'`

## If NovaSeq is contained in the description field of the sample sheet then set the pixel distance appropriately
## Assumption: all read groups come from some sequencer model. otherwise pixel distance would be set to NovaSeq
## If mixing NovaSeq and non-NovaSeq then this workflow would need to change.

if [[ $SEQUENCER_MODEL == *"NovaSeq"* ]]
	then
		PIXEL_DISTANCE="2500"
	else
		PIXEL_DISTANCE="100"
fi

## --Mark Duplicates with Picard, write a duplicate report

START_MARK_DUPLICATES=`date '+%s'`

	# if any part of pipe fails set exit to non-zero

	set -o pipefail

	$JAVA_1_8/java -jar \
		-Xmx16g \
		-XX:ParallelGCThreads=4 \
		$PICARD_DIR/picard.jar \
		MarkDuplicates \
		ASSUME_SORT_ORDER=queryname \
		$INPUT \
		OUTPUT=/dev/stdout \
		VALIDATION_STRINGENCY=SILENT \
		METRICS_FILE=$CORE_PATH/$PROJECT/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" \
		COMPRESSION_LEVEL=0 \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=$PIXEL_DISTANCE \
	| $SAMBAMBA_DIR/sambamba \
		sort \
		-t 4 \
		-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".dup.bam" \
		/dev/stdin

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

END_MARK_DUPLICATES=`date '+%s'`

echo $SM_TAG"_"$PROJECT",C.01,MARK_DUPLICATES,"$HOSTNAME","$START_MARK_DUPLICATES","$END_MARK_DUPLICATES \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar \
-Xmx16g \
-XX:ParallelGCThreads=4 \
$PICARD_DIR/picard.jar \
MarkDuplicates \
ASSUME_SORT_ORDER=queryname \
$INPUT \
OUTPUT=/dev/stdout \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=$CORE_PATH/$PROJECT/REPORTS/PICARD_DUPLICATES/$SM_TAG"_MARK_DUPLICATES.txt" \
COMPRESSION_LEVEL=0 \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=$PIXEL_DISTANCE \
\| $SAMBAMBA_DIR/sambamba \
sort \
-t 4 \
-o $CORE_PATH/$PROJECT/TEMP/$SM_TAG".dup.bam" \
/dev/stdin \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

# exit with the signal from the program

	exit $SCRIPT_STATUS
