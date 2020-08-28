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

	BAIT_BED=$5

		BAIT_BED_NAME=(`basename $BAIT_BED .bed`)

	TARGET_BED=$6

		TARGET_BED_NAME=(`basename $TARGET_BED .bed`)

	TITV_BED=$7

		TITV_BED_NAME=(`basename $TITV_BED_NAME .bed`)

	REF_GENOME=$8
		REF_DIR=$(dirname $REF_GENOME)
		REF_BASENAME=$(basename $REF_GENOME | sed 's/.fasta//g ; s/.fa//g')

# FIX BED FILES (FOR GRCH37)

	# FIX THE BAIT BED FILE

		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# remove CHR PREFIXES (THIS IS FOR GRCH37)
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
					
			awk 1 $BAIT_BED | sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"BAIT_BED_NAME".bed"

	# FIX THE TARGET BED FILE

		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# remove CHR PREFIXES (THIS IS FOR GRCH37)
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
					
			awk 1 $TARGET_BED | sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"TARGET_BED_NAME".bed"

	# FIX THE TITV BED FILE

		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# remove CHR PREFIXES (THIS IS FOR GRCH37)
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
					
			awk 1 $TITV_BED | sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"TITV_BED_NAME".bed"

# MAKE PICARD INTERVAL FILES (1-based start)

	# bait bed

		(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
			; awk 1 $BAIT_BED \
				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
				| awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}') \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnBait.picard.bed"

	# target bed

		(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
			; awk 1 $TARGET_BED \
				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
				| awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}') \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTarget.picard.bed"
