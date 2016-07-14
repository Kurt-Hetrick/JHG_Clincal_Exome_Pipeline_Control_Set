# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash


# tell sge to submit any of these queue when available
#$ -q rnd.q,prod.q,test.q

# tell sge that you are in the users current working directory
#$ -cwd

# tell sge to export the users environment variables
#$ -V

# tell sge to submit at this priority setting
#$ -p -10

# tell sge to output both stderr and stdout to the same file
#$ -j y

# export all variables, useful to find out what compute node the program was executed on
# redirecting stderr/stdout to file as a log.

set

JAVA_1_7=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
FAMILY=$5
REF_GENOME=$6


# Filter to just on all of the variants all

START_FILTER_TO_FAMILY_VARIANTS_PASS=`date '+%s'`

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--keepOriginalAC \
--sample_file $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".sample.list" \
--variant $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.vcf" \
--excludeNonVariants \
--excludeFiltered \
-o $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.VARIANT_ONLY.PASS.vcf"

END_FILTER_TO_FAMILY_VARIANTS_PASS=`date '+%s'`

HOSTNAME=`hostname`

echo $FAMILY"_"$PROJECT",N.001,FILTER_TO_FAMILY_VARIANTS_PASS,"$HOSTNAME","$START_FILTER_TO_FAMILY_VARIANTS_PASS","$END_FILTER_TO_FAMILY_VARIANTS_PASS \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--keepOriginalAC \
--sample_file $CORE_PATH/$PROJECT/$FAMILY/$FAMILY/$FAMILY".sample.list" \
--variant $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.vcf" \
--excludeNonVariants \
--excludeFiltered \
-o $CORE_PATH/$PROJECT/$FAMILY/VCF/VQSR/$FAMILY".VQSR.ANNOTATED.JUST_FAMILY.VARIANT_ONLY.PASS.vcf" \
>> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/$FAMILY/$FAMILY".COMMAND.LINES.txt"
