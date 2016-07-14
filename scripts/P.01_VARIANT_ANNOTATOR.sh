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

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3
PED_FILE=$4

PROJECT=$5
REF_GENOME=$6


START_ADD_MORE_ANNOTATION=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--annotation AlleleBalance \
--annotation AlleleBalanceBySample \
--annotation AlleleCountBySample \
--annotation GCContent \
--annotation GenotypeSummaries \
--annotation HardyWeinberg \
--annotation HomopolymerRun \
--annotation InbreedingCoeff \
--annotation MVLikelihoodRatio \
--annotation SampleList \
--annotation TandemRepeatAnnotator \
--annotation TransmissionDisequilibriumTest \
--annotation VariantType \
--resource:OneKGP /isilon/sequencing/1000genomes/Full_Project/Sep_2014/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf \
--expression OneKGP.AF \
--expression OneKGP.EAS_AF \
--expression OneKGP.AMR_AF \
--expression OneKGP.AFR_AF \
--expression OneKGP.EUR_AF \
--expression OneKGP.SAS_AF \
--resourceAlleleConcordance \
--pedigree $PED_FILE \
--pedigreeValidationType SILENT \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.VQSR.vcf \
-L $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.VQSR.vcf \
-o $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz

END_ADD_MORE_ANNOTATION=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",M.001,ADD_MORE_ANNOTATION,"$HOSTNAME","$START_ADD_MORE_ANNOTATION","$END_ADD_MORE_ANNOTATION \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $REF_GENOME \
--disable_auto_index_creation_and_locking_when_reading_rods \
--annotation AlleleBalance \
--annotation AlleleBalanceBySample \
--annotation AlleleCountBySample \
--annotation GCContent \
--annotation GenotypeSummaries \
--annotation HardyWeinberg \
--annotation HomopolymerRun \
--annotation InbreedingCoeff \
--annotation MVLikelihoodRatio \
--annotation SampleList \
--annotation TandemRepeatAnnotator \
--annotation TransmissionDisequilibriumTest \
--annotation VariantType \
--resource:OneKGP /isilon/sequencing/1000genomes/Full_Project/Sep_2014/20130502/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf \
--expression OneKGP.AF \
--expression OneKGP.EAS_AF \
--expression OneKGP.AMR_AF \
--expression OneKGP.AFR_AF \
--expression OneKGP.EUR_AF \
--expression OneKGP.SAS_AF \
--resourceAlleleConcordance \
--pedigree $PED_FILE \
--pedigreeValidationType SILENT \
--variant $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.VQSR.vcf \
-L $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.VQSR.vcf \
-o $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz \
>> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

md5sum $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

md5sum $CORE_PATH/$PROJECT/JOINT_VCF/CONTROL_DATA_SET.VQSR.ANNOTATED.vcf.gz.tbi \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
