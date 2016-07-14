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

PROJECT=$4
REF_GENOME=$5
DBSNP=$6

RIS_ID=${SM_TAG%@*}
BARCODE_2D=${SM_TAG#*@}

CONTROL_REPO="/isilon/sequencing/Seq_Proj/CGC_CONTROL_DATA_SET"

START_GENOTYPE_GVCF=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--standard_min_confidence_threshold_for_calling 30 \
--standard_min_confidence_threshold_for_emitting 0 \
--annotation AS_BaseQualityRankSumTest \
--annotation AS_FisherStrand \
--annotation AS_InbreedingCoeff \
--annotation AS_MappingQualityRankSumTest \
--annotation AS_RMSMappingQuality \
--annotation AS_ReadPosRankSumTest \
--annotation AS_StrandOddsRatio \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--annotation LikelihoodRankSumTest \
--excludeIntervals 1:145017822-145017822 \
--variant $CONTROL_REPO/Control_GVCF.list \
-o $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf

END_GENOTYPE_GVCF=`date '+%s'`

HOSTNAME=`hostname`

echo $PROJECT",I.001,GENOTYPE_GVCF,"$HOSTNAME","$START_GENOTYPE_GVCF","$END_GENOTYPE_GVCF \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R $REF_GENOME \
--dbsnp $DBSNP \
--annotateNDA \
--includeNonVariantSites \
--disable_auto_index_creation_and_locking_when_reading_rods \
--standard_min_confidence_threshold_for_calling 30 \
--standard_min_confidence_threshold_for_emitting 0 \
--annotation AS_BaseQualityRankSumTest \
--annotation AS_FisherStrand \
--annotation AS_InbreedingCoeff \
--annotation AS_MappingQualityRankSumTest \
--annotation AS_RMSMappingQuality \
--annotation AS_ReadPosRankSumTest \
--annotation AS_StrandOddsRatio \
--annotation FractionInformativeReads \
--annotation StrandBiasBySample \
--annotation StrandAlleleCountsBySample \
--annotation LikelihoodRankSumTest \
--excludeIntervals 1:145017822-145017822 \
--variant $CONTROL_REPO/Control_GVCF.list \
-o $CORE_PATH/$PROJECT/TEMP/CONTROL_DATA_SET.RAW.vcf \
>> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt

echo >> $CORE_PATH/$PROJECT/CONTROL_DATA_SET.COMMAND.LINES.txt
