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
	SM_TAG=$4

	CODING_BED=$5
		CODING_BED_NAME=(`basename $CODING_BED .bed`)
	TARGET_BED=$6
		TARGET_BED_NAME=(`basename $TARGET_BED .bed`)
	BAIT_BED=$7
		BAIT_BED_NAME=(`basename $BAIT_BED .bed`)
	TITV_BED=$8
		TITV_BED_NAME=(`basename $TITV_BED_NAME .bed`)


# FIX THE TARGET BED FILE
	# make sure that there is EOF
	# remove CARRIAGE RETURNS
	# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
# PAD THE REFSEQ CODING BED FILE BY 10 BASES. NEED TO CHECK WHAT THIS IS USED FOR LATER.

	awk 1 $CODING_BED \
		| sed 's/\r//g' \
		| sed -r 's/[[:space:]]+/\t/g' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10}' \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_CODING.bed"

# FIX THE TARGET BED FILE
	# make sure that there is EOF
	# remove CARRIAGE RETURNS
	# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
# PAD THE TARGET BED FILE BY X BP (i THINK THIS WILL BE DEFINED BY GUI)
	# THIS IS FOR SLICING

	awk 1 $TARGET_BED \
		| sed 's/\r//g' \
		| sed -r 's/[[:space:]]+/\t/g' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10}' \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_TARGET.bed"

# FIX THE BAIT BED FILE
	# make sure that there is EOF
	# remove CARRIAGE RETURNS
	# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
# FOR DATA PROCESSING AND METRICS REPORTS

	awk 1 $BAIT_BED \
		| sed 's/\r//g' \
		| sed -r 's/[[:space:]]+/\t/g' \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed"

# PAD THE BAIT FILE.
# THE BAIT FILE IS THE COMBINED MERGING OF THE CIDR TWIST BAIT BED FILE
## AND THE CODING BED FILE WHICH IS REFSEQ SELECT CDS AND MISSING OMIM.
# THIS WILL BE USED FOR GVCF FILE CREATION SO IT WILL BE SUPER PADDED.
# MERGE THE PADDED THE TARGET BED WITH THE BAIT BED FILE

	cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_TARGET.bed" \
	$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"BAIT_BED_NAME".bed" \
		| sort -k 1,1 -k 2,2n -k 3,3n \
		| singularity exec $ALIGNMENT_CONTAINER bedtools merge -i - \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_GVCF.bed"

######################################################################################

# 	REF_GENOME=$8
# 		REF_DIR=$(dirname $REF_GENOME)
# 		REF_BASENAME=$(basename $REF_GENOME | sed 's/.fasta//g ; s/.fa//g')

# # MAKE PICARD INTERVAL FILES (1-based start)

# 	# bait bed

# 		(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
# 			; awk 1 $BAIT_BED \
# 				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
# 				| awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}') \
# 		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnBait.picard.bed"

# 	# target bed

# 		(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
# 			; awk 1 $TARGET_BED \
# 				| sed -r 's/\r//g ; s/chr//g ; s/[[:space:]]+/\t/g' \
# 				| awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}') \
# 		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".OnTarget.picard.bed"
