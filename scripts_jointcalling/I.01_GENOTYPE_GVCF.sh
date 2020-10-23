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
	CONTROL_REPO=$6

START_GENOTYPE_GVCF=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
--logging_level ERROR \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--variant $CONTROL_REPO/Control_GVCF.list \
-o $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf

# --excludeIntervals 1:145017822-145017822 \
# removing this for the moment b/c I want to see how the June 26th nightly handles it.

END_GENOTYPE_GVCF=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",I.001,GENOTYPE_GVCF,"$HOSTNAME","$START_GENOTYPE_GVCF","$END_GENOTYPE_GVCF \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
--logging_level ERROR \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--variant $CONTROL_REPO/Control_GVCF.list \
-o $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf \
>> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt
