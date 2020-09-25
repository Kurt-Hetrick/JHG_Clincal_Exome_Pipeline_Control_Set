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
		CODING_BED_NAME=$(basename $CODING_BED .bed)
		CODING_MD5=$(md5sum $CODING_BED | cut -c 1-7)
	# note: since the coding bed file is supposed to be static.
	# i'm tracking it by creating a md5 and applying it to the file name
	# which gets captured by the command line output
	# i don't care how the bed file is tracked when changes are made, just that it is tracked.
	TARGET_BED=$6
		TARGET_BED_NAME=$(basename $TARGET_BED .bed)
	BAIT_BED=$7
		BAIT_BED_NAME=$(basename $BAIT_BED .bed)
	TITV_BED=$8
		TITV_BED_NAME=$(basename $TITV_BED .bed)
	CYTOBAND_BED=$9
	REF_GENOME=${10}
		REF_DIR=$(dirname $REF_GENOME)
		REF_BASENAME=$(basename $REF_GENOME | sed 's/.fasta//g ; s/.fa//g')
	PADDING_LENGTH=${11}

# FIX AND PAD THE CODING BED FILE
	# make sure that there is EOF
	# remove CARRIAGE RETURNS
	# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
	# PAD THE REFSEQ CODING BED FILE BY THE PADDING LENGTH
	# remove chr prefix
	# remove MT genome (done in another pipeline)
	# remove annotation fields

		awk 1 $CODING_BED \
			| sed 's/\r//g' \
			| sed -r 's/[[:space:]]+/\t/g' \
			| awk 'BEGIN {OFS="\t"} {print $1,$2-"'$PADDING_LENGTH'",$3+"'$PADDING_LENGTH'"}' \
			| sed 's/^chr//g' \
			| grep -v "^MT" \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$CODING_BED_NAME"-"$CODING_MD5"-"$PADDING_LENGTH"-BP-PAD.bed"

# FIX AND PAD THE TARGET BED FILE
	# make sure that there is EOF
	# remove CARRIAGE RETURNS
	# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
	# PAD THE TARGET BED FILE BY THE PADDING LENGTH
	# remove chr prefix
	# remove MT genome (done in another pipeline)
	# THIS IS FOR SLICING

		awk 1 $TARGET_BED \
			| sed 's/\r//g' \
			| sed -r 's/[[:space:]]+/\t/g' \
			| awk 'BEGIN {OFS="\t"} {print $1,$2-"'$PADDING_LENGTH'",$3+"'$PADDING_LENGTH'"}' \
			| sed 's/^chr//g' \
			| grep -v "^MT" \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$TARGET_BED_NAME"-"$PADDING_LENGTH"-BP-PAD.bed"

# FIX THE CODING BED FILE. THIS IS TO BE COMBINED WITH THE CODING BED FILE
# AND THEN PADDED BY 250 BP AND THEN MERGED (FOR OVERLAPPING INTERVALS) FOR GVCF CREATION.
	# FOR DATA PROCESSING AND METRICS REPORTS AS WELL.
		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
		# remove chr prefix
		# remove MT genome (done in another pipeline)

			awk 1 $CODING_BED \
				| sed 's/\r//g' \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/^chr//g' \
				| grep -v "^MT" \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$CODING_BED_NAME"-"$CODING_MD5".bed"

# FIX AND PAD THE ANNOTATED CODING BED FILE
	# make sure that there is EOF
	# remove CARRIAGE RETURNS
	# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
	# remove chr prefix
	# remove MT genome (done in another pipeline)

		awk 1 $CODING_BED \
			| sed 's/\r//g' \
			| sed -r 's/[[:space:]]+/\t/g' \
			| sed 's/^chr//g' \
			| grep -v "^MT" \
			| awk 'BEGIN {OFS="\t"} {print $1,$2-"'$PADDING_LENGTH'",$3+"'$PADDING_LENGTH'",$4,$5,$6,$7}' \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_"$CODING_BED_NAME"-"$CODING_MD5"-"$PADDING_LENGTH"-BP-PAD-ANNOTATED.bed"

# FIX THE BAIT BED FILE. THIS IS TO BE COMBINED WITH THE BAIT BED FILE AND THEN PADDED BY 250 BP
# AND THEN MERGED (FOR OVERLAPPING INTERVALS) FOR GVCF CREATION.
	# FOR DATA PROCESSING AND METRICS REPORTS AS WELL.
		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
		# remove chr prefix
		# remove MT genome (done in another pipeline)

			awk 1 $BAIT_BED \
				| sed 's/\r//g' \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/^chr//g' \
				| grep -v "^MT" \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed"

# FIX THE TITV BED FILE FOR DATA PROCESSING AND METRICS REPORTS.
		# make sure that there is EOF
		# remove CARRIAGE RETURNS
		# CONVERT VARIABLE LENGTH WHITESPACE FIELD DELIMETERS TO SINGLE TAB.
		# remove chr prefix
		# remove MT genome (done in another pipeline)

			awk 1 $TITV_BED \
				| sed 's/\r//g' \
				| sed -r 's/[[:space:]]+/\t/g' \
				| sed 's/^chr//g' \
				| grep -v "^MT" \
			>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$TITV_BED_NAME".bed"

# THE GVCF BED FILE IS THE COMBINED MERGING OF THE CIDR TWIST BAIT BED FILE
## AND THE CODING BED FILE WHICH IS REFSEQ SELECT CDS AND MISSING OMIM.
# THIS WILL BE USED FOR GVCF FILE CREATION SO IT WILL BE SUPER PADDED WITH 250 BP.

	cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$CODING_BED_NAME"-"$CODING_MD5".bed" \
	$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed" \
		| sort -k 1,1 -k 2,2n -k 3,3n \
		| awk 'BEGIN {OFS="\t"} {print $1,$2-"'$PADDING_LENGTH'",$3+"'$PADDING_LENGTH'"}' \
		| singularity exec $ALIGNMENT_CONTAINER bedtools merge -i - \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_GVCF.bed"

# Format the cytoband file.
# strip out the "chr" prefix from the chromsome name
# print the chromsome, start, end, the first character of the cytoband (to get the chromosome arm).
# the file is already sorted correctly so group by chromosome and chromosome arm and print the first start and last end
	# for the chromosome/arm combination
# print CHROMOSOME, START, END, ARM (TAB DELIMITED) TO MAKE A BED FILE.

	sed 's/^chr//g' $CYTOBAND_BED \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,substr($4,0,1)}' \
		| singularity exec $ALIGNMENT_CONTAINER datamash \
			-s \
			-g 1,4 \
			first 2 \
			last 3 \
		| awk 'BEGIN {OFS="\t"} {print $1,$3,$4,$2}' \
	>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".CHROM_ARM.bed"

# MAKE PICARD INTERVAL FILES (1-based start)
# ti/tv bed is used as the target since it shouldn't change
	# GRAB THE SEQUENCING DICTIONARY FORM THE ".dict" file in the directory where the reference genome is located
	# then concatenate with the fixed bed file.
	# add 1 to the start
	# picard interval needs strand information and a locus name
		# made everything plus stranded b/c i don't think this information is used
		# constructed locus name with chr name, start+1, stop

	# bait bed

		(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
			; awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' \
				$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME".bed") \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$BAIT_BED_NAME"-picard.bed"

	# target-TITV bed

		(grep "^@SQ" $REF_DIR/$REF_BASENAME".dict" \
			; awk 'BEGIN {OFS="\t"} {print $1,($2+1),$3,"+",$1"_"($2+1)"_"$3}' \
				$CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$TITV_BED_NAME".bed") \
		>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"$TITV_BED_NAME"-picard.bed"
