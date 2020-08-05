# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash


# tell sge to submit any of these queue when available
#$ -q bigdata.q,c6320.q,lemon.q,prod.q,rnd.q,c6420_21.q,c6420_23.q

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

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3

PROJECT=$4
SM_TAG=$5
REF_GENOME=$6
TARGET_BED=$7


# Filter to just on all of the variants all

START_FILTER_TO_SAMPLE_VARIANT_ONLY_PASS_TARGET=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--keepOriginalAC \
-L $TARGET_BED \
--sample_name $SM_TAG \
--excludeNonVariants \
--excludeFiltered \
--keepOriginalDP \
--removeUnusedAlternates \
--variant $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz \
-o $CORE_PATH/$PROJECT/VCF/FILTERED_ON_TARGET/$SM_TAG".VARIANT_SITES.PASS.ON_TARGET.vcf.gz"

END_FILTER_TO_SAMPLE_VARIANT_ONLY_PASS_TARGET=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",S.001,FILTER_TO_SAMPLE_VARIANT_ONLY_PASS_TARGET,"$HOSTNAME","$START_FILTER_TO_SAMPLE_VARIANT_ONLY_PASS_TARGET","$END_FILTER_TO_SAMPLE_VARIANT_ONLY_PASS_TARGET \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R $REF_GENOME \
--keepOriginalAC \
-L $TARGET_BED \
--sample_name $SM_TAG \
--excludeNonVariants \
--excludeFiltered \
--keepOriginalDP \
--removeUnusedAlternates \
--variant $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz \
-o $CORE_PATH/$PROJECT/VCF/FILTERED_ON_TARGET/$SM_TAG".VARIANT_SITES.PASS.ON_TARGET.vcf.gz" \
>> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

md5sum $CORE_PATH/$PROJECT/VCF/FILTERED_ON_TARGET/$SM_TAG".VARIANT_SITES.PASS.ON_TARGET.vcf.gz" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/VCF/FILTERED_ON_TARGET/$SM_TAG".VARIANT_SITES.PASS.ON_TARGET.vcf.gz.tbi" \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

