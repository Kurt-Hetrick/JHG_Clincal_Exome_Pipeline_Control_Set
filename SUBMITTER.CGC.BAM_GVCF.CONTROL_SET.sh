#!/bin/bash

SAMPLE_SHEET=$1
PED_FILE=$2

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

SCRIPT_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINES/JHGenomics_CGC_Clinical_Exome_Control_Set/scripts"
# The above hash value is the corresponding commit at https://github.com/Kurt-Hetrick/JHGenomics_CGC_Clinical_Exome_Control_Set

CORE_PATH="/mnt/clinical/ddl/NGS/Exome_Data"
CONTROL_REPO="/mnt/clinical/ddl/NGS/Exome_Resources/CONTROL_REPO"

# PIPELINE PROGRAMS
JAVA_1_6="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/jre1.6.0_25/bin"
JAVA_1_8="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/jdk1.8.0_73/bin"
BWA_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/bwa-0.7.8"
PICARD_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/picard-tools-2.1.1"
GATK_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/GenomeAnalysisTK-3.7"
VERIFY_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/verifyBamID_20120620/bin/"
TABIX_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/tabix-0.2.6"
SAMTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/samtools-0.1.18"
DATAMASH_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/datamash-1.0.6"
BEDTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/bedtools-2.22.0/bin"
VCFTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/vcftools_0.1.12b/bin"
PLINK2_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/PLINK2"
KING_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/KING/Linux-king19"
CIDRSEQSUITE_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/CIDRSeqSuiteSoftware_Version_4_0/"
ANNOVAR_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/ANNOVAR/2013_09_11"

# PIPELINE FILES
GENE_LIST="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/RefSeqGene.GRCh37.Ready.txt"
VERIFY_VCF="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
CODING_BED="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/RefSeq.Unique.GRCh37.FINAL.19Feb2018.bed"
CYTOBAND_BED="/isilon/cgc/PIPELINE_FILES/GRCh37.Cytobands.bed"
HAPMAP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/hapmap_3.3.b37.vcf"
OMNI_1KG="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_omni2.5.b37.vcf"
HI_CONF_1KG_PHASE1_SNP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_phase1.snps.high_confidence.b37.vcf"
MILLS_1KG_GOLD_INDEL="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf"
PHASE3_1KG_AUTOSOMES="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
DBSNP_129="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/dbsnp_138.b37.excluding_sites_after_129.vcf"

##### MAKE A DIRECTORY TREE ##### SHOULD BE COMPLETE #####

mkdir -p ~/CGC_PIPELINE_TEMP

MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
PED_PREFIX=`basename $PED_FILE .ped`

##########################################################

# I typically comment out the bed file making after I make them the first time.

SETUP_PROJECT ()
{
FORMAT_MANIFEST
MERGE_PED_MANIFEST
CREATE_SAMPLE_INFO_ARRAY
MAKE_PROJ_DIR_TREE
echo "echo Making padded annotated RefSeq coding bed file for $SAMPLE"
PAD_REFSEQ
# echo "echo Making padded target bed file for $SAMPLE"
# PAD_TARGET
# echo "echo Making everything merged together bait file for $SAMPLE"
# MAKE_BAIT
}

FORMAT_MANIFEST ()
{
sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| sed 's/,/\t/g' \
| sort -k 8,8 \
>| ~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt
}

MERGE_PED_MANIFEST ()
{
awk 1 $PED_FILE \
| sed 's/\r//g' \
| sort -k 2,2 \
| join -1 8 -2 2 -e '-'  -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.1,2.3,2.4,2.5,2.6' \
~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt /dev/stdin \
>| ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt
}

# MAKE AN ARRAY FOR EACH SAMPLE
	## SAMPLE_INFO_ARRAY[0] = PROJECT
	## SAMPLE_INFO_ARRAY[1] = FAMILY
	## SAMPLE_INFO_ARRAY[2] = SM_TAG
		## SAMPLE = SM_TAG
	## SAMPLE_INFO_ARRAY[3] = BAIT BED FILE
	## SAMPLE_INFO_ARRAY[4] = TARGET_BED_FILE

CREATE_SAMPLE_INFO_ARRAY ()
{
SAMPLE_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$20,$8,$15,$16}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
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
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/DEPTH_OF_COVERAGE/{TARGET,REFSEQ_CODING_PLUS_10bp} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/INSERT_SIZE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/LOCAL_REALIGNMENT_INTERVALS \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/REPORTS/ANEUPLOIDY_CHECK \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/JOINT_VCF/ \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}_ANNOVAR \
$CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/{FASTQ,LOGS,COMMAND_LINES}
}

# PAD THE REFSEQ canonical transcript bed file by 10 bases.
# can make this as an input variable with a default value 10 if i have to ever give more than 0 effs.

PAD_REFSEQ ()
{
awk 1 $CODING_BED \
| sed 's/\r//g' \
| sed -r 's/[[:space:]]+/\t/g' \
| awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10}' \
>| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}"_PADDED_CODING.bed"
}

# PAD THE TARGET BED FILE BY 10 BP

PAD_TARGET ()
{
awk 1 ${SAMPLE_INFO_ARRAY[4]} \
| sed 's/\r//g' \
| sed -r 's/[[:space:]]+/\t/g' \
| awk 'BEGIN {OFS="\t"} {print $1,$2-10,$3+10}' \
>| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}"_PADDED_TARGET.bed"
}

# MERGE THE PADDED THE TARGET BED WITH THE BAIT BED FILE

MAKE_BAIT ()
{
cat $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}"_PADDED_TARGET.bed" \
${SAMPLE_INFO_ARRAY[3]} \
| sort -k 1,1 -k 2,2n -k 3,3n \
| $BEDTOOLS_DIR/bedtools merge -i - \
>| $CORE_PATH/${SAMPLE_INFO_ARRAY[0]}/TEMP/${SAMPLE_INFO_ARRAY[2]}_BAIT.bed
}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
SETUP_PROJECT
done

############################################################

# to create the qsub cmd line to submit bwa alignments to the cluster
# handle blank lines
# handle something else too

awk 'BEGIN {FS="\t"} {split($19,INDEL,";");split($8,smtag,"[@-]"); \
print "qsub","-N","A.01_BWA_"$8"_"$2"_"$3"_"$4,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$8"_"$2"_"$3"_"$4".BWA.log",\
"'$SCRIPT_DIR'""/A.01_BWA.sh",\
"'$BWA_DIR'","'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12"\n""sleep 1s"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt

# create a hold job id qsub command line based on the number of
# submit merging the bam files created by bwa mem above
# only launch when every lane for a sample is done being processed by bwa mem

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam"}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1,2 collapse 3 collapse 4 \
| awk 'BEGIN {FS="\t"} \
gsub(/,/,",A.01_BWA_"$2"_",$3) \
gsub(/,/,",INPUT=" "'$CORE_PATH'" "/" $1 "/TEMP/",$4) \
{print "qsub","-N","B.01_MERGE_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".MERGE.BAM.FILES.log",\
"-hold_jid","A.01_BWA_"$2"_"$3, \
"'$SCRIPT_DIR'""/B.01_MERGE_SORT_AGGRO.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'",$1,$2,"INPUT=" "'$CORE_PATH'" "/" $1 "/TEMP/" $4 "\n" "sleep 1s"}'

# Mark duplicates on the bam file above. Create a Mark Duplicates report which goes into the QC report

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$19}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$19}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$19,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","F.01_PERFORM_BQSR_"$2"_"$1,\
"-hold_jid","E.01_INDEL_REALIGNER_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PERFORM_BQSR.log",\
"'$SCRIPT_DIR'""/F.01_PERFORM_BQSR.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,INDEL[1],INDEL[2],$5"\n""sleep 1s"}'

# write Final Bam file

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.01_HAPLOTYPE_CALLER_"$1"_"$2,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".HAPLOTYPE_CALLER.log",\
"'$SCRIPT_DIR'""/H.01_HAPLOTYPE_CALLER.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# Run POST BQSR TABLE

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$19,$18}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
print "qsub","-N","H.02-A.01_ANALYZE_COVARIATES_"$2"_"$1,\
"-hold_jid","H.02_POST_BQSR_TABLE_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".ANALYZE_COVARIATES.log",\
"'$SCRIPT_DIR'""/H.02-A.01_ANALYZE_COVARIATES.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# RUN DOC RefSeq CODING PLUS 10 BP FLANKS
# This will specifically be for the RefSeg transcript IDs merged with UCSC canonical transcripts

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".DOC_CODING_10bpFLANKS.log",\
"'$SCRIPT_DIR'""/H.03_DOC_CODING_10bpFLANKS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3"\n""sleep 1s"}'

# RUN ANEUPLOIDY_CHECK AFTER DOC RefSeq CODING PLUS 10 BP FLANKS
# redo

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.03-A.01_DOC_CHROM_DEPTH_"$2"_"$1,\
"-hold_jid","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".ANEUPLOIDY_CHECK.log",\
"'$SCRIPT_DIR'""/H.03-A.01_CHROM_DEPTH.sh",\
"'$CORE_PATH'","'$CYTOBAND_BED'","'$DATAMASH_DIR'","'$BEDTOOLS_DIR'",$1,$2"\n""sleep 1s"}'

# RUN FORMATTING PER BASE COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{split($2,smtag,"[@]"); \
print "qsub", "-N", "H.03-A.02_PER_BASE_" smtag[1] "_" smtag[2] "_" $1,\
"-hold_jid", "H.03_DOC_CODING_10bpFLANKS_" $2 "_" $1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE.log",\
"'$SCRIPT_DIR'""/H.03-A.02_PER_BASE.sh",\
"'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2"\n""sleep 1s"}'

# RUN FILTERING PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{split($2,smtag,"[@]"); \
print "qsub","-N","H.03-A.02-A.01_PER_BASE_FILTER_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE_FILTER.log",\
"'$SCRIPT_DIR'""/H.03-A.02-A.01_PER_BASE_FILTERED.sh",\
"'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{split($2,smtag,"[@]"); \
print "qsub","-N","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE_BGZIP.log",\
"'$SCRIPT_DIR'""/H.03-A.02-A.02_PER_BASE_BGZIP.sh",\
"'$CORE_PATH'","'$TABIX_DIR'",$1,$2"\n""sleep 1s"}'

# TABIX PER BASE COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{split($2,smtag,"[@]"); \
print "qsub","-N","H.03-A.02-A.02-A.01_PER_BASE_TABIX_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE_TABIX.log",\
"'$SCRIPT_DIR'""/H.03-A.02-A.02-A.01_PER_BASE_TABIX.sh",\
"'$CORE_PATH'","'$TABIX_DIR'",$1,$2"\n""sleep 1s"}'

# RUN FORMATTING PER CODING INTERVAL COVERAGE WITH GENE NAME ANNNOTATION

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{split($2,smtag,"[@]"); \
print "qsub","-N","H.03-A.03_PER_INTERVAL_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_INTERVAL.log",\
"'$SCRIPT_DIR'""/H.03-A.03_PER_INTERVAL.sh",\
"'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2"\n""sleep 1s"}'

# RUN FILTERING PER CODING INTERVAL COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{split($2,smtag,"[@]"); \
print "qsub","-N","H.03-A.03_PER_INTERVAL_FILTER_"smtag[1]"_"smtag[2]"_"$1,\
"-hold_jid","H.03-A.03_PER_INTERVAL_"smtag[1]"_"smtag[2]"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_INTERVAL_FILTER.log",\
"'$SCRIPT_DIR'""/H.03-A.03-A.01_PER_INTERVAL_FILTERED.sh",\
"'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# RUN DOC TARGET BED (Generally this with all RefGene coding exons unless it becomes targeted)

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.05_DOC_TARGET_BED_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".DOC_TARGET_BED.log",\
"'$SCRIPT_DIR'""/H.05_DOC_TARGET_BED.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3"\n""sleep 1s"}'

# RUN COLLECT MULTIPLE METRICS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$18,$15}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{split($2,smtag,"[@-]"); \
print "qsub","-N","H.07_COLLECT_HS_METRICS_"$2"_"$1,\
"-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".COLLECT_HS_METRICS.log",\
"'$SCRIPT_DIR'""/H.07_COLLECT_HS_METRICS.sh",\
"'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3"\n""sleep 1s"}'

# RUN SELECT VERIFYBAM ID VCF

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$15}' \
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

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
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
SAMPLE_INFO_ARRAY_VERIFY_BAM=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$20,$8,$12,$15}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
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
#     1  CGC_CONTROL_SET_3_7_REFSEQ_TEMP
#     2  HV3JGBCXX
#     3  1
#     4  CTGAGCCA
#     5  ILLUMINA
#     6  D01_CRE10-1_H05
#     7  7/29/2016
#     8  CRE10-1
#     9  Johns_Hopkins_DNA_Diagnostic_Lab
#    10  HiSeq2500_RapidRun
#    11  HV3JGBCXX_1_CTGAGCCA_D01_CRE10-1_H05
#    12  /mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/human_g1k_v37_decoy.fasta
#    13  MBS
#    14  -2
#    15  /isilon/cgc/BED_FILES/RefGene.Coding.Clean.Format.bed
#    16  /isilon/cgc/BED_FILES/ALLBED_BED_File_Agilent_ClinicalExome_S06588914_RefSeqCoding_051017_noCHR.bed
#    17  /isilon/cgc/BED_FILES/RefGene.Coding.Clean.Format.bed
#    18  /mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/dbsnp_138.b37.vcf
#    19  /mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_phase1.indels.b37.vcf;/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf
#    20  CRE10
#    21  0
#    22  0
#    23  1
#    24  2
#
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
