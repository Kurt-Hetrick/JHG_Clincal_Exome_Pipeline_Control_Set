# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash


# tell sge to submit any of these queue when available
#$ -q bigdata.q,lemon.q,prod.q,rnd.q,uhoh.q

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
REF_GENOME=$5

START_APPLY_RECALIBRATION_SNP=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf \
--ts_filter_level 99.9 \
-recalFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.recal \
-tranchesFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.tranches \
-mode SNP \
-o $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.VQSR.SNP.vcf

END_APPLY_RECALIBRATION_SNP=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",K.001,APPLY_RECALIBRATION_SNP,"$HOSTNAME","$START_APPLY_RECALIBRATION_SNP","$END_APPLY_RECALIBRATION_SNP \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--input:VCF $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf \
--ts_filter_level 99.9 \
-recalFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.recal \
-tranchesFile $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.HC.SNV.tranches \
-mode SNP \
-o $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.VQSR.SNP.vcf \
>> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt
