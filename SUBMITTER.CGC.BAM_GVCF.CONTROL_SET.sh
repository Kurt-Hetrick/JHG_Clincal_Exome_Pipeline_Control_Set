#!/bin/bash

module load sge

SAMPLE_SHEET=$1
PED_FILE=$2

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

SCRIPT_DIR="/isilon/sequencing/Kurt/GIT_REPO/JHGenomics_CGC_Clinical_Exome_Control_Set/scripts"

JAVA_1_8="/isilon/sequencing/Kurt/Programs/Java/jdk1.8.0_73/bin"
CORE_PATH="/isilon/sequencing/Seq_Proj/"
BWA_DIR="/isilon/sequencing/Kurt/Programs/BWA/bwa-0.7.8"
PICARD_DIR="/isilon/sequencing/Kurt/Programs/Picard/picard-tools-2.1.1"
GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-nightly-2016-07-12-gaa9ac69"
VERIFY_DIR="/isilon/sequencing/Kurt/Programs/VerifyBamID/verifyBamID_20120620/bin/"
GENE_LIST="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/RefSeqGene.GRCh37.Ready.txt"
VERIFY_VCF="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
CODING_BED="/isilon/sequencing/CIDRSeqSuiteSoftware/RELEASES/5.0.0/aux_files/UCSC_hg19_CodingOnly_083013_MERGED_noContigs_noCHR.bed"
SAMTOOLS_DIR="/isilon/sequencing/Kurt/Programs/samtools/samtools-0.1.18"
TABIX_DIR="/isilon/sequencing/Kurt/Programs/TABIX/tabix-0.2.6"
CORE_PATH="/isilon/sequencing/Seq_Proj"
DATAMASH_DIR="/isilon/sequencing/Kurt/Programs/PATH"
CYTOBAND_BED="/isilon/sequencing/Kurt/CGC/GRCh37.Cytobands.bed"
# BEDTOOLS IS v2.22.0
BEDTOOLS_DIR="/isilon/sequencing/Kurt/Programs/PATH"

##### MAKE A DIRECTORY TREE ##### SHOULD BE COMPLETE #####

mkdir -p ~/CGC_PIPELINE_TEMP

MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
PED_PREFIX=`basename $PED_FILE .ped`

##########################################################

SETUP_PROJECT ()
{
FORMAT_MANIFEST
MERGE_PED_MANIFEST
CREATE_SAMPLE_INFO_ARRAY
MAKE_PROJ_DIR_TREE
}

FORMAT_MANIFEST ()
{
sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| sed 's/,/\t/g' \
| sort -k 8 \
>| ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt
}

MERGE_PED_MANIFEST ()
{
awk 1 $PED_FILE \
| sed 's/\r//g' \
| sort -k 2 \
| join -1 8 -2 2 ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt /dev/stdin \
| awk 'BEGIN {OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$1,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' \
>| ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt
}

# MAKE AN ARRAY FOR EACH SAMPLE
	## SAMPLE_INFO_ARRAY[0] = PROJECT
	## SAMPLE_INFO_ARRAY[1] = FAMILY
	## SAMPLE_INFO_ARRAY[2] = SM_TAG
		## SAMPLE = SM_TAG

CREATE_SAMPLE_INFO_ARRAY ()
{
SAMPLE_INFO_ARRAY=(`awk '$8=="'$SAMPLE'" {print $1,$19,$8}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

# PROJECT DIRECTORY TREE CREATOR

MAKE_PROJ_DIR_TREE ()
{
mkdir -p $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/BAM \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/HC_BAM \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/INDEL/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/SNV/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/MIXED/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/GVCF \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,VERIFYBAMID_CHR} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/CONCORDANCE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/COUNT_COVARIATES/{GATK_REPORT,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/DEPTH_OF_COVERAGE/{TARGET,UCSC_CODING_PLUS_10bp} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/INSERT_SIZE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/ANEUPLOIDY_CHECK \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/JOINT_VCF/ \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/{TEMP,FASTQ,LOGS,COMMAND_LINES}
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
SETUP_PROJECT
done

############################################################

# to create the qsub cmd line to submit bwa alignments to the cluster
# handle blank lines
# handle something else too

awk '{split($18,INDEL,";");split($8,smtag,"[@-]"); \
print "qsub","-N","A.01_BWA_"$8"_"$2"_"$3"_"$4,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$8"_"$2"_"$3"_"$4".BWA.log",\
"'$SCRIPT_DIR'""/A.01_BWA.sh",\
"'$BWA_DIR'","'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12"\n""sleep 1s"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt

# create a hold job id qsub command line based on the number of
# submit merging the bam files created by bwa mem above
# only launch when every lane for a sample is done being processed by bwa mem

awk 'BEGIN {OFS="\t"} {print $1,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1,2 collapse 3 collapse 4 \
| awk 'BEGIN {FS="\t"} \
gsub(/,/,",A.01_BWA_"$2"_",$3) \
gsub(/,/,",INPUT=/isilon/sequencing/Seq_Proj/"$1"/TEMP/",$4) \
{print "qsub","-N","B.01_MERGE_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".MERGE.BAM.FILES.log",\
"-hold_jid","A.01_BWA_"$2"_"$3, \
"'$SCRIPT_DIR'""/B.01_MERGE_SORT_AGGRO.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,"INPUT=/isilon/sequencing/Seq_Proj/"$1"/TEMP/"$4"\n""sleep 1s"}'

# Mark duplicates on the bam file above. Create a Mark Duplicates report which goes into the QC report

awk 'BEGIN {OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","C.01_MARK_DUPLICATES_"$2"_"$1,\
"-hold_jid","B.01_MERGE_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".MARK_DUPLICATES.log",\
"'$SCRIPT_DIR'""/C.01_MARK_DUPLICATES.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# Generate a list of places that could be potentially realigned.

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","D.01_REALIGNER_TARGET_CREATOR_"$2"_"$1,\
"-hold_jid","C.01_MARK_DUPLICATES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".REALIGNER_TARGET_CREATOR.log",\
"'$SCRIPT_DIR'""/D.01_REALIGNER_TARGET_CREATOR.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,INDEL[1],INDEL[2]"\n""sleep 1s"}'

# With the list generated above walk through the BAM file and realign where necessary
# Write out a new bam file

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","E.01_INDEL_REALIGNER_"$2"_"$1,\
"-hold_jid","D.01_REALIGNER_TARGET_CREATOR_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".INDEL_REALIGNER.log",\
"'$SCRIPT_DIR'""/E.01_INDEL_REALIGNER.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,INDEL[1],INDEL[2]"\n""sleep 1s"}'

# Run Base Quality Score Recalibration

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$18,$17,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","F.01_PERFORM_BQSR_"$2"_"$1,\
"-hold_jid","E.01_INDEL_REALIGNER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PERFORM_BQSR.log",\
"'$SCRIPT_DIR'""/F.01_PERFORM_BQSR.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,INDEL[1],INDEL[2],$5,$6"\n""sleep 1s"}'

# write Final Bam file

awk 'BEGIN {OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","G.01_FINAL_BAM_"$2"_"$1,\
"-hold_jid","F.01_PERFORM_BQSR_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FINAL_BAM.log",\
"'$SCRIPT_DIR'""/G.01_FINAL_BAM.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

##### ALL H.00X SERIES OF SCRIPTS CAN BE RUN IN PARALLEL SINCE THEY ARE DEPENDENT ON FINAL BAM FILE GENERATION #####

# Run Haplotype Caller in GVCF mode

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.01_HAPLOTYPE_CALLER_"$1"_"$2,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".HAPLOTYPE_CALLER.log",\
"'$SCRIPT_DIR'""/H.01_HAPLOTYPE_CALLER.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 1s"}'

# Run POST BQSR TABLE

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$18,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","H.02_POST_BQSR_TABLE_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".POST_BQSR_TABLE.log",\
"'$SCRIPT_DIR'""/H.02_POST_BQSR_TABLE.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,INDEL[1],INDEL[2],$5"\n""sleep 1s"}'

# Run ANALYZE COVARIATES

awk 'BEGIN {OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","H.02-A.01_ANALYZE_COVARIATES_"$2"_"$1,\
"-hold_jid","H.02_POST_BQSR_TABLE_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".ANALYZE_COVARIATES.log",\
"'$SCRIPT_DIR'""/H.02-A.01_ANALYZE_COVARIATES.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# RUN DOC CODING PLUS 10 BP FLANKS

awk 'BEGIN {OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".DOC_CODING_10bpFLANKS.log",\
"'$SCRIPT_DIR'""/H.03_DOC_CODING_10bpFLANKS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$CODING_BED'","'$GENE_LIST'",$1,$2,$3"\n""sleep 1s"}'

# RUN ANEUPLOIDY_CHECK AFTER DOC TARGET BED FINISHES

awk 'BEGIN {OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.03-A.01_DOC_CHROM_DEPTH_"$2"_"$1,\
"-hold_jid","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".ANEUPLOIDY_CHECK.log",\
"'$SCRIPT_DIR'""/H.03-A.01_CHROM_DEPTH.sh",\
"'$CORE_PATH'","'$CYTOBAND_BED'","'$DATAMASH_DIR'","'$BEDTOOLS_DIR'",$1,$2"\n""sleep 1s"}'

# RUN DOC TARGET BED

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.05_DOC_TARGET_BED_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".DOC_TARGET_BED.log",\
"'$SCRIPT_DIR'""/H.05_DOC_TARGET_BED.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3,$4"\n""sleep 1s"}'

# RUN COLLECT MULTIPLE METRICS

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$17,$14}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.06_COLLECT_MULTIPLE_METRICS_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".COLLECT_MULTIPLE_METRICS.log",\
"'$SCRIPT_DIR'""/H.06_COLLECT_MULTIPLE_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# RUN COLLECT HS METRICS

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$15,$16}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.07_COLLECT_HS_METRICS_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".COLLECT_HS_METRICS.log",\
"'$SCRIPT_DIR'""/H.07_COLLECT_HS_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# RUN SELECT VERIFYBAM ID VCF

awk 'BEGIN {OFS="\t"} {print $1,$8,$12,$14}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.08_SELECT_VERIFYBAMID_VCF_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".SELECT_VERIFYBAMID_VCF.log",\
"'$SCRIPT_DIR'""/H.08_SELECT_VERIFYBAMID_VCF.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$VERIFY_VCF'",$1,$2,$3,$4"\n""sleep 1s"}'

# RUN VERIFYBAMID ALL

awk 'BEGIN {OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.08-A.01_VERIFYBAMID_"$2"_"$1,\
"-hold_jid","H.08_SELECT_VERIFYBAMID_VCF_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".VERIFYBAMID.log",\
"'$SCRIPT_DIR'""/H.08-A.01_VERIFYBAMID.sh",\
"'$CORE_PATH'","'$VERIFY_DIR'",$1,$2"\n""sleep 1s"}'

###################################################
### RUN VERIFYBAM ID PER CHROMOSOME - VITO ########
###################################################

CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM ()
{
SAMPLE_INFO_ARRAY_VERIFY_BAM=(`awk '$8=="'$SAMPLE'" {print $1,$19,$8,$12,$14}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
}

CALL_SELECT_VERIFY_BAM ()
{
echo \
qsub \
-N H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
-hold_jid G.01_FINAL_BAM_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.SELECT_VERIFYBAMID_chr$CHROMOSOME.log \
$SCRIPT_DIR/H.09_SELECT_VERIFYBAMID_VCF_CHR.sh \
$JAVA_1_8 $GATK_DIR $CORE_PATH $VERIFY_VCF \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[3]} \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[4]} $CHROMOSOME
}

CALL_VERIFYBAMID ()
{
echo \
qsub \
-N H.09-A.01_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
-hold_jid H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.VERIFYBAMID_chr$CHROMOSOME.log \
$SCRIPT_DIR/H.09-A.01_VERIFYBAMID_CHR.sh \
$CORE_PATH $VERIFY_DIR \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} \
$CHROMOSOME
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
	for CHROMOSOME in {1..22}
		do
		CALL_SELECT_VERIFY_BAM
		echo sleep 1s
		CALL_VERIFYBAMID
		echo sleep 1s
	done
done

#####################################################
### JOIN THE PER CHROMOSOME VERIFYBAMID REPORTS #####
#####################################################

BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR ()
{
	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
	do
	HOLD_ID_PATH="-hold_jid "
	for CHROMOSOME in {{1..22},{X,Y}};
 	do
 		HOLD_ID_PATH=$HOLD_ID_PATH"H.09-A.01_VERIFYBAMID_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}"_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}"_"chr$CHROMOSOME","
 	done
 done
}

 CAT_VERIFYBAMID_CHR ()
 {
echo \
qsub \
-N H.09-A.01-A.01_JOIN_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
$HOLD_ID_PATH \
-o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.CAT_VERIFYBAMID_CHR.log \
$SCRIPT_DIR/H.09-A.01-A.01_CAT_VERIFYBAMID_CHR.sh \
$CORE_PATH \
${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}
 }

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
 do
 	CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
	BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR
	CAT_VERIFYBAMID_CHR
	echo sleep 1s
 done

#############################################


### kEY FOR BLAH ###
#
#      1  CGC_160212_HJLWVBCXX_CGCDev1_TEST
#      2  HJLWVBCXX
#      3  1
#      4  ATGCCTAA
#      5  ILLUMINA
#      6  A01_NA12878
#      7  2/12/2016
#      8  NA12878
#      9  CGC
#     10  HiSeq2500_RapidRun
#     11  HJLWVBCXX_1_ATGCCTAA_A01_NA12878
#     12  /isilon/sequencing/GATK_resource_bundle/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.fasta
#     13  MBS
#     14  /isilon/sequencing/data/Work/BED/Production_BED_files/TsTv_BED_File_Agilent_ClinicalExome_S06588914_OnExon_merged_021015_noCHR.bed
#     15  /isilon/sequencing/data/Work/BED/Production_BED_files/ALLBED_BED_File_Agilent_ClinicalExome_S06588914_ALLBed_merged_021015_noCHR.bed
#     16  /isilon/sequencing/data/Work/BED/Production_BED_files/Targets_BED_File_Agilent_ClinicalExome_S06588914_OnTarget_merged_noCHR_013015.bed
#     17  /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.vcf
#     18  /isilon/sequencing/GATK_resource_bundle/2.2/b37/1000G_phase1.indels.b37.vcf;/isilon/sequencing/GATK_resource_bundle/2.2/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
#     19  XC01463
#     20  NA12891
#     21  NA12892
#     22  2
#     23  2
#######


###### SAMPLE MANIFEST KEY...NOT SURE WHAT I AM GOING TO END UP DOING HERE ######

# PROJECT=$1 # the Seq Proj folder name. 1st column in sample manifest
# FLOWCELL=$2 # flowcell that sample read group was performed on. 2nd column of sample manifest
# LANE=$3 # lane of flowcell that sample read group was performed on. 3rd column of the sample manifest
# INDEX=$4 # sample barcode. 4th column of the sample manifest
# PLATFORM=$5 # type of sequencing chemistry matching SAM specification. 5th column of the sample manifest.
# LIBRARY_NAME=$6 # library group of the sample read group.
# 								# Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not
# 								# 6th column of the sample manifest
# RUN_DATE=$7 # should be the run set up date to match the seq run folder name, but it has been arbitrarily populated. field X of manifest.
# SM_TAG=$8 # sample ID. sample name for all files, etc. field X of manifest
# CENTER=$9 # the center/funding mechanism. field X of manifest.
# DESCRIPTION=${10} # Generally we use to denote the sequencer setting (e.g. rapid run). field X of manifest.
# REF_GENOME=${11} # the reference genome used in the analysis pipeline. field X of manifest.
# TI_TV_BED=${12} # populated from sample manifest. where ucsc coding exons overlap with bait and target bed files
# BAIT_BED=${13} # populated from sample manifest. a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
# 								# Used for limited where to run base quality score recalibration on where to create gvcf files.
# TARGET_BED=${14} # populated from sample manifest. bed file acquired from manufacturer of their targets. field X of sample manifest.
# DBSNP=${15} # populated from sample manifest. used to annotate ID field in VCF file. masking in base call quality score recalibration.
# KNOWN_INDEL_1=${16} # populated from sample manifest. used for BQSR masking, sensitivity in local realignment.
# KNOWN_INDEL_2=${17} # populated from sample manifest. used for BQSR masking, sensitivity in local realignment.
#
# RIS_ID=${SM_TAG%@*} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
# BARCODE_2D=${SM_TAG#*@} # no longer needed when using PHOENIX. used to needed to break out the "@" in the sm tag so it wouldn't break things.
#
####################################################################################

#### BOILERPLATE...I HAVE NOT DECIDED WHAT I AM GOING TO DO HERE######

# function GRAB_MANIFEST {
# sed 's/\r//g' $SAMPLE_SHEET \
# | awk 'BEGIN {FS=","} NR>1 \
# {split($19,INDEL,";");split($8,smtag,"@");print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2]}'
# }
#
# function GRAB_PROJECT_NAMES {
# PROJECT_NAMES=`sed 's/\r//g' $SAMPLE_SHEET \
# | awk 'BEGIN {FS=","} NR>1 print $1}'`
# }
######################################################################
