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

	ANNOVAR_DIR=$1
	CORE_PATH=$2

	PROJECT=$3
	SM_TAG=$4

# if any part of pipe fails set exit to non-zero

	set -o pipefail

# Grab the header lines staring with #
# Grab the column descriptor line and "FILTER" to the end.
# Run annovar to extract the entire vcf row for each allele and add it to the end of the canonical first 5 rows.
# Do a first five column join b/w the original annovar report and the .annovar just created, appending the FILTER field to the end of the original ANNOVAR report

	grep "^#" $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"/$SM_TAG".VARIANT_SITES"/$SM_TAG".VARIANT_SITES_ANNOVAR_REPORT.txt" \
	>| $CORE_PATH/$PROJECT/REPORTS/ANNOVAR/$SM_TAG".VARIANT_SITES_ANNOVAR_REPORT.txt" ; \
	grep -v "^#" $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"/$SM_TAG".VARIANT_SITES"/$SM_TAG".VARIANT_SITES_ANNOVAR_REPORT.txt" \
		| awk 'NR==1 {print $0"\t""FILTER"}' \
	>> $CORE_PATH/$PROJECT/REPORTS/ANNOVAR/$SM_TAG".VARIANT_SITES_ANNOVAR_REPORT.txt" ; \
	$ANNOVAR_DIR/convert2annovar.pl \
		--format vcf4 \
	$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"/$SM_TAG".VARIANT_SITES.vcf" \
		--includeinfo \
		| awk 'NR==FNR{a[$1,$2,$3,$4,$5]=$12;next} ($1,$2,$3,$4,$5) in a{print $0 "\t" a[$1,$2,$3,$4,$5]}' \
		/dev/stdin \
		$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"/$SM_TAG".VARIANT_SITES"/$SM_TAG".VARIANT_SITES_ANNOVAR_REPORT.txt" \
	>> $CORE_PATH/$PROJECT/REPORTS/ANNOVAR/$SM_TAG".VARIANT_SITES_ANNOVAR_REPORT.txt"

# check the exit signal at this point.

	SCRIPT_STATUS=`echo $?`

# delete annovar staging area

	# rm -rvf $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"/

# exit with the signal from the program

	exit $SCRIPT_STATUS
