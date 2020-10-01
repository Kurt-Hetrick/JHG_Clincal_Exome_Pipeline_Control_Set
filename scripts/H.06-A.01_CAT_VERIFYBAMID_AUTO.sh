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

	CORE_PATH=$1
	DATAMASH_DIR=$2

	PROJECT=$3
	SM_TAG=$4
	TARGET_BED=$5
		TARGET_BED_NAME=(`basename $TARGET_BED .bed`)

echo \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt"

for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG"-"TARGET_BED_NAME".bed" \
	| sed -r 's/[[:space:]]+/\t/g' \
	| cut -f 1 \
	|  egrep -v "X|Y|MT" \
	| sort \
	| uniq \
	| $DATAMASH_DIR/datamash collapse 1 \
	| sed 's/,/ /g');
	do
		# I'm stripping out the "chr" prefix here b/c I don't want to deal with it...I should specify the column in case SM_TAG contain chr...
		cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG"."$CHROMOSOME".selfSM" \
			| grep -v ^# \
			| awk 'BEGIN {OFS="\t"} {print($1,"'$CHROMOSOME'",$7,$4,$8,$9,$6)}' \
			| sed 's/chr//g' \
		>> $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt"
done

sed -i '/^\s*$/d' $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt"

(awk '$2~/^[0-9]/' $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt" | sort -k2,2n ; \
awk '$2=="X"' $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt" ; \
awk '$2=="Y"' $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt" ; \
awk '$2=="MT"' $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_unsorted.txt") \
>| $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_joined.txt"

echo "#SM_TAG" CHROM VERIFYBAM_FREEMIX VERIFYBAM_SNPS VERIFYBAM_FREELK1 VERRIFYBAM_FREELK0 VERIFYBAM_AVG_DP \
>| $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID_CHR/$SM_TAG".VERIFYBAMID.PER_CHR.txt"

cat $CORE_PATH/$PROJECT/TEMP/$SM_TAG".verifybamID_joined.txt" \
>> $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID_CHR/$SM_TAG".VERIFYBAMID.PER_CHR.txt"

sed -i 's/ /\t/g' $CORE_PATH/$PROJECT/REPORTS/VERIFYBAMID_CHR/$SM_TAG".VERIFYBAMID.PER_CHR.txt"
