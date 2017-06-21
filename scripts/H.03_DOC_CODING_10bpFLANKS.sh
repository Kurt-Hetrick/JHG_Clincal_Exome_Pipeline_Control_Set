# ---qsub parameter settings---
# --these can be overrode at qsub invocation--

# tell sge to execute in bash
#$ -S /bin/bash

# tell sge to submit any of these queue when available
#$ -q cgc.q

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

echo

JAVA_1_8=$1
GATK_DIR=$2
CORE_PATH=$3
# CODING_BED=$4
GENE_LIST=$4

PROJECT=$5
SM_TAG=$6
REF_GENOME=$7

### --Depth of Coverage ALL REFSEQ CODING EXONS PLUS 10 BP FLANKS--
### the padding should be an input variable
### Ben will probably be creating a padded bed file in the very beginning.
### actually i'm doing the same thing, just not using it here.

START_DOC_CODING=`date '+%s'`

$JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $REF_GENOME \
-I $CORE_PATH/$PROJECT/BAM/$SM_TAG".bam" \
-geneList:REFSEQ $GENE_LIST \
-L $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_CODING.bed" \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-o $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks" \
-ct 10 \
-ct 15 \
-ct 20 \
-ct 30 \
-ct 50

END_DOC_CODING=`date '+%s'`

HOSTNAME=`hostname`

echo $SM_TAG"_"$PROJECT",H.001,DOC_CODING,"$HOSTNAME","$START_DOC_CODING","$END_DOC_CODING \
>> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".WALL.CLOCK.TIMES.csv"

echo $JAVA_1_8/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $REF_GENOME \
-I $CORE_PATH/$PROJECT/BAM/$SM_TAG".bam" \
-geneList:REFSEQ $GENE_LIST \
-L $CORE_PATH/$PROJECT/TEMP/$SM_TAG"_PADDED_CODING.bed" \
-mmq 20 \
-mbq 10 \
--outputFormat csv \
-o $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks" \
-ct 10 \
-ct 15 \
-ct 20 \
-ct 30 \
-ct 50 \
>> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

echo >> $CORE_PATH/$PROJECT/COMMAND_LINES/$SM_TAG".COMMAND.LINES.txt"

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.EveryBase.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.EveryBase.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_cumulative_coverage_counts" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_cumulative_coverage_counts.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_cumulative_coverage_counts.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_cumulative_coverage_proportions" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_cumulative_coverage_proportions.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_cumulative_coverage_proportions.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_gene_summary" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_gene_summary.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_gene_summary.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_statistics" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_statistics.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_statistics.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_summary" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_summary.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_interval_summary.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_statistics" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_statistics.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_statistics.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"

#####

mv -v $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_summary" \
$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_summary.csv"

# md5sum $CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/REFSEQ_CODING_PLUS_10bp/$SM_TAG".ALL_REFSEQ_CODING_10bpFlanks.sample_summary.csv" \
# >> $CORE_PATH/$PROJECT/REPORTS/$PROJECT".CIDR.Analysis.MD5.txt"
