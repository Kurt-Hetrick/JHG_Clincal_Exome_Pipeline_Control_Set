#!/usr/bin/env bash

# INPUT VARIABLES

	SAMPLE_SHEET=$1
	PED_FILE=$2
	PADDING_LENGTH=$3 # optional. if no 3rd argument present then the default is 10
	# THIS PAD IS FOR SLICING

		if [[ ! $PADDING_LENGTH ]]
			then
			PADDING_LENGTH="10"
		fi

	QUEUE_LIST=$4 # optional. if no 4th argument present then the default is cgc.q
		# if you want to set this then you need to set the 3rd argument as well (even to the default)

		if [[ ! $QUEUE_LIST ]]
			then
			QUEUE_LIST="cgc.q"
		fi

	PRIORITY=$5 # optional. if no 5th argument present then the default is -15.
		# if you want to set this then you need to set the 3rd and 4th argument as well (even to the default)

			if [[ ! $PRIORITY ]]
				then
				PRIORITY="-15"
			fi

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

	SUBMITTER_SCRIPT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

	SCRIPT_DIR="$SUBMITTER_SCRIPT_PATH/scripts_jointcalling"

##################
# CORE VARIABLES #
##################

	## This will always put the current working directory in front of any directory for PATH
	## added /bin for RHEL6

		export PATH=".:$PATH:/bin"

	# where the input/output sequencing data will be located.

		CORE_PATH="/mnt/clinical/ddl/NGS/Exome_Data"

	# WHERE THE CONTROL DATA SET RESIDES

		CONTROL_REPO="/mnt/clinical/ddl/NGS/Exome_Resources/TWIST_CONTROL_SET1.200601_PIPELINE_2_0_0"

	# used for tracking in the read group header of the cram file

		PIPELINE_VERSION=`git --git-dir=$SCRIPT_DIR/../.git --work-tree=$SCRIPT_DIR/.. log --pretty=format:'%h' -n 1`
	# bind the host file system /mnt to the singularity container. in case I use it in the submitter.

		export SINGULARITY_BINDPATH="/mnt:/mnt"

	# QSUB ARGUMENTS LIST
		# set shell on compute node
		# start in current working directory
		# transfer submit node env to compute node
		# set SINGULARITY BINDPATH
		# set queues to submit to
		# set priority
		# combine stdout and stderr logging to same output file

			QSUB_ARGS="-S /bin/bash" \
				QSUB_ARGS=$QSUB_ARGS" -cwd" \
				QSUB_ARGS=$QSUB_ARGS" -V" \
				QSUB_ARGS=$QSUB_ARGS" -v SINGULARITY_BINDPATH=/mnt:/mnt" \
				QSUB_ARGS=$QSUB_ARGS" -q $QUEUE_LIST" \
				QSUB_ARGS=$QSUB_ARGS" -p $PRIORITY" \
				QSUB_ARGS=$QSUB_ARGS" -j y"

#####################
# PIPELINE PROGRAMS #
#####################

	ALIGNMENT_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/images/ddl_ce_control_align-0.0.3.simg"
	# contains the following software and is on Ubuntu 16.04.5 LTS
		# gatk 4.0.11.0 (base image). also contains the following.
			# Python 3.6.2 :: Continuum Analytics, Inc.
				# samtools 0.1.19
				# bcftools 0.1.19
				# bedtools v2.25.0
				# bgzip 1.2.1
				# tabix 1.2.1
				# samtools, bcftools, bgzip and tabix will be replaced with newer versions.
				# R 3.2.5
					# dependencies = c("gplots","digest", "gtable", "MASS", "plyr", "reshape2", "scales", "tibble", "lazyeval")    # for ggplot2
					# getopt_1.20.0.tar.gz
					# optparse_1.3.2.tar.gz
					# data.table_1.10.4-2.tar.gz
					# gsalib_2.1.tar.gz
					# ggplot2_2.2.1.tar.gz
				# openjdk version "1.8.0_181"
				# /gatk/gatk.jar -> /gatk/gatk-package-4.0.11.0-local.jar
		# added
			# picard.jar 2.17.0 (as /gatk/picard.jar)
			# samblaster-v.0.1.24
			# sambamba-0.6.8
			# bwa-0.7.15
			# datamash-1.6
			# verifyBamID v1.1.3
			# samtools 1.10
			# bgzip 1.10
			# tabix 1.10
			# bcftools 1.10.2

	GATK_3_7_0_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/images/gatk3-3.7-0.simg"
	# singularity pull docker://broadinstitute/gatk3:3.7-0
	# used for generating the depth of coverage reports.
		# comes with R 3.1.1 with appropriate packages needed to create gatk pdf output
		# also comes with some version of java 1.8
		# jar file is /usr/GenomeAnalysisTK.jar

	# PIPELINE PROGRAMS

		SAMTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/samtools-0.1.18"
		VCFTOOLS_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/vcftools_0.1.12b/bin"
		PLINK2_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/PLINK2"
		KING_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/KING/Linux-king19"

		JAVA_1_6="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/jre1.6.0_25/bin"
		CIDRSEQSUITE_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/CIDRSeqSuiteSoftware_Version_4_0/"
		ANNOVAR_DIR="/mnt/clinical/ddl/NGS/Exome_Resources/PROGRAMS/ANNOVAR/2013_09_11"

##################
# PIPELINE FILES #
##################

	GENE_LIST="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/RefSeqGene.GRCh37.rCRS.MT.bed"
		# md5 dec069c279625cfb110c2e4c5480e036
	VERIFY_VCF="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
	CODING_BED="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINES/TWIST/JHGenomics_CGC_Clinical_Exome_Control_Set/GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_HGNC_annotated.bed"
	CYTOBAND_BED="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/GRCh37.Cytobands.bed"
	HAPMAP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/hapmap_3.3.b37.vcf"
	OMNI_1KG="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_omni2.5.b37.vcf"
	HI_CONF_1KG_PHASE1_SNP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_phase1.snps.high_confidence.b37.vcf"
	MILLS_1KG_GOLD_INDEL="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf"
	PHASE3_1KG_AUTOSOMES="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
	DBSNP_129="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/dbsnp_138.b37.excluding_sites_after_129.vcf"

#################################
##### MAKE A DIRECTORY TREE #####
#################################

	mkdir -p ~/CGC_PIPELINE_TEMP

	MANIFEST_PREFIX=`basename $SAMPLE_SHEET .csv`
	PED_PREFIX=`basename $PED_FILE .ped`

	FORMAT_MANIFEST ()
	{
		awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
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
		| join -1 8 -2 2 -e '-'  -t $'\t' \
		-o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.1,2.3,2.4,2.5,2.6' \
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

	CREATE_SAMPLE_ARRAY ()
	{
			SAMPLE_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" \
				{split($19,INDEL,";"); \
				print $1,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],$20,$21,$22,$23,$24}' \
					~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
					| sort \
					| uniq`)

			#  1  Project=the Seq Proj folder name

				PROJECT=${SAMPLE_ARRAY[0]}

			################################################################################
			# 2 SKIP : FCID=flowcell that sample read group was performed on ###############
			# 3 SKIP : Lane=lane of flowcell that sample read group was performed on] ######
			# 4 SKIP : Index=sample barcode ################################################
			# 5 SKIP : Platform=type of sequencing chemistry matching SAM specification ####
			# 6 SKIP : Library_Name=library group of the sample read group #################
			# 7 SKIP : Date=should be the run set up date to match the seq run folder name #
			################################################################################

			#  8  SM_Tag=sample ID

				SM_TAG=${SAMPLE_ARRAY[1]}
					SGE_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g') # If there is an @ in the qsub or holdId name it breaks

			#  9  Center=the center/funding mechanism

				CENTER=${SAMPLE_ARRAY[2]}

			# 10  Description=Generally we use to denote the sequencer setting (e.g. rapid run)
			# “HiSeq-X”, “HiSeq-4000”, “HiSeq-2500”, “HiSeq-2000”, “NextSeq-500”, or “MiSeq”.

				SEQUENCER_MODEL=${SAMPLE_ARRAY[3]}

			#########################
			# 11  SKIP : Seq_Exp_ID #
			#########################

			# 12  Genome_Ref=the reference genome used in the analysis pipeline

				REF_GENOME=${SAMPLE_ARRAY[4]}

			#####################################
			# 13  Operator: SKIP ################
			# 14  Extra_VCF_Filter_Params: SKIP #
			#####################################

			# 15  TS_TV_BED_File=where ucsc coding exons overlap with bait and target bed files

				TITV_BED=${SAMPLE_ARRAY[5]}

			# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
			# Used for limited where to run base quality score recalibration on where to create gvcf files.

				BAIT_BED=${SAMPLE_ARRAY[6]}

			# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.

				TARGET_BED=${SAMPLE_ARRAY[7]}

			# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in base call quality score recalibration.

				DBSNP=${SAMPLE_ARRAY[8]}

			# 19  KNOWN_INDEL_FILES=used for BQSR masking, sensitivity in local realignment.

				KNOWN_INDEL_1=${SAMPLE_ARRAY[9]}
				KNOWN_INDEL_2=${SAMPLE_ARRAY[10]}

			# 20 family that sample belongs to

				FAMILY=${SAMPLE_ARRAY[11]}

			# 21 MOM

				MOM=${SAMPLE_ARRAY[12]}

			# 22 DAD

				DAD=${SAMPLE_ARRAY[13]}

			# 23 GENDER

				GENDER=${SAMPLE_ARRAY[14]}

			# 24 PHENOTYPE

				PHENOTYPE=${SAMPLE_ARRAY[15]}
	}

# PROJECT DIRECTORY TREE CREATOR

	MAKE_PROJ_DIR_TREE ()
	{
		mkdir -p $CORE_PATH/$PROJECT/{FASTQ,LOGS,COMMAND_LINES,CRAM,HC_CRAM,GVCF,JOINT_VCF} \
		$CORE_PATH/$PROJECT/INDEL/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/SNV/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/MIXED/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/VCF/{FILTERED_ON_BAIT,FILTERED_ON_TARGET} \
		$CORE_PATH/$PROJECT/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,VERIFYBAMID_AUTO,ANEUPLOIDY_CHECK,RG_HEADER,QUALITY_YIELD,ERROR_SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT \
		$CORE_PATH/$PROJECT/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/{TARGET_PADDED,CODING_PADDED} \
		$CORE_PATH/$PROJECT/REPORTS/HYB_SELECTION/PER_TARGET_COVERAGE \
		$CORE_PATH/$PROJECT/REPORTS/INSERT_SIZE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/REPORTS/MEAN_QUALITY_BY_CYCLE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/TEMP/$SM_TAG"_ANNOVAR"
	}

	SETUP_PROJECT ()
	{
		FORMAT_MANIFEST
		MERGE_PED_MANIFEST
		CREATE_SAMPLE_ARRAY
		MAKE_PROJ_DIR_TREE
	}

for SAMPLE in $(awk 1 $SAMPLE_SHEET \
		| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
		| awk 'BEGIN {FS=","} NR>1 {print $8}' \
		| sort \
		| uniq );
	do
		SETUP_PROJECT
done

################################
#### JOINT CALLING AND VQSR ####
################################

### Run GenotypeGVCF per Family

	GENOTYPE_GVCF_ALL_CONTROLS ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N I.01_GENOTYPE_GVCF"_"$PROJECT \
			-o $CORE_PATH/$PROJECT/LOGS/$PROJECT".GENOTYPE_GVCF.log" \
		$SCRIPT_DIR/I.01_GENOTYPE_GVCF.sh \
			$GATK_3_7_0_CONTAINER \
			$CORE_PATH \
			$PROJECT \
			$SM_TAG \
			$BAIT_BED \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

### Run Variant Recalibrator for the SNP model, this is done in parallel with the INDEL model

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12,$18}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","J.01_VARIANT_RECALIBRATOR_SNP_"$1,\
"-hold_jid","I.01_GENOTYPE_GVCF_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".VARIANT_RECALIBRATOR_SNP.log",\
"'$SCRIPT_DIR'""/J.01_VARIANT_RECALIBRATOR_SNP.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,"'$HAPMAP'","'$OMNI_1KG'","'$HI_CONF_1KG_PHASE1_SNP'""\n""sleep 3s"}'

### Run Variant Recalibrator for the INDEL model, this is done in parallel with the SNP model

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","J.02_VARIANT_RECALIBRATOR_INDEL_"$1,\
"-hold_jid","I.01_GENOTYPE_GVCF_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".VARIANT_RECALIBRATOR_INDEL.log",\
"'$SCRIPT_DIR'""/J.02_VARIANT_RECALIBRATOR_INDEL.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,"'$MILLS_1KG_GOLD_INDEL'""\n""sleep 3s"}'

### Run Apply Recalbration with the SNP model to the VCF file

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","K.01_APPLY_RECALIBRATION_SNP_"$1,\
"-hold_jid","J.01_VARIANT_RECALIBRATOR_SNP_"$1",""J.02_VARIANT_RECALIBRATOR_INDEL_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".APPLY_RECALIBRATION_SNP.log",\
"'$SCRIPT_DIR'""/K.01_APPLY_RECALIBRATION_SNP.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

### Run Apply Recalibration with the INDEL model to the VCF file.

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","L.01_APPLY_RECALIBRATION_INDEL_"$1,\
"-hold_jid","K.01_APPLY_RECALIBRATION_SNP_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".APPLY_RECALIBRATION_INDEL.log",\
"'$SCRIPT_DIR'""/L.01_APPLY_RECALIBRATION_INDEL.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

### Add all possible GATK annotations to the VCF file.

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","P.01_VARIANT_ANNOTATOR_"$1,\
"-hold_jid","L.01_APPLY_RECALIBRATION_INDEL_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".VARIANT_ANNOTATOR.log",\
"'$SCRIPT_DIR'""/P.01_VARIANT_ANNOTATOR.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$PED_FILE'",$1,$2,"'$PHASE3_1KG_AUTOSOMES'""\n""sleep 3s"}'

##### DOING VCF BREAKOUTS #####

### SUBSETTING FROM COHORT (FAMILY PLUS CONTROL SET) VCF ###

# FILTER TO JUST VARIANT SITES

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","S.01_FILTER_COHORT_VARIANT_ONLY_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".FILTER_COHORT_VARIANT_ONLY.log",\
"'$SCRIPT_DIR'""/S.01_FILTER_COHORT_VARIANT_ONLY.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

# FILTER TO JUST PASSING VARIANT SITES

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","S.02_FILTER_COHORT_VARIANT_ONLY_PASS_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".FILTER_COHORT_VARIANT_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.02_FILTER_COHORT_VARIANT_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

### SUBSETTING TO SAMPLE VCFS ###

## SUBSET TO SAMPLE VCF ALL SITES ##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.06_FILTER_TO_SAMPLE_ALL_SITES_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.06_FILTER_TO_SAMPLE_ALL_SITES.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSET TO SAMPLE VARIANTS ONLY 

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.07_FILTER_TO_SAMPLE_VARIANTS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_VARIANTS.log",\
"'$SCRIPT_DIR'""/S.07_FILTER_TO_SAMPLE_VARIANTS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING VARIANTS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.08_FILTER_TO_SAMPLE_VARIANTS_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_VARIANTS_PASS.log",\
"'$SCRIPT_DIR'""/S.08_FILTER_TO_SAMPLE_VARIANTS_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING SNVS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09_FILTER_TO_SNV_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_SNV_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.09_FILTER_TO_SAMPLE_SNV_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING INDELS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.10_FILTER_TO_INDEL_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_INDEL_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.10_FILTER_TO_SAMPLE_INDEL_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSET TO SAMPLE PASSING MIXED

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.11_FILTER_TO_MIXED_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_MIXED_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.11_FILTER_TO_SAMPLE_MIXED_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 3s"}'

## SUBSET TO TARGET SNV ONLY PASS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_TARGET_SNV_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.12_FILTER_TO_SAMPLE_TARGET_SNV_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO TARGET INDEL ONLY PASS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_TARGET_INDEL_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.13_FILTER_TO_SAMPLE_TARGET_INDEL_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO TARGET MIXED ONLY PASS

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_TARGET_MIXED_ONLY_PASS.log",\
"'$SCRIPT_DIR'""/S.14_FILTER_TO_SAMPLE_TARGET_MIXED_ONLY_PASS.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

## SUBSET TO SAMPLE VCF ALL SITES ON TARGET##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$17}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET_"$2"_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_ALL_SITES_TARGET.log",\
"'$SCRIPT_DIR'""/S.15_FILTER_TO_SAMPLE_ALL_SITES_TARGET.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'


########################
##### TITV SECTION #####
########################

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"$2"_"$1,\
"-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_TITV_VCF.log",\
"'$SCRIPT_DIR'""/S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4"\n""sleep 3s"}'

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND OVERLAP WITH DBSNP 129

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"$2"_"$1,\
"-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_TITV_VCF_KNOWN.log",\
"'$SCRIPT_DIR'""/S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,"'$DBSNP_129'""\n""sleep 3s"}'

# BREAK DOWN TO ALL PASSING SNV THAT FALL IN TITV BED FILE AND DO NOT OVERLAP WITH DBSNP 129

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$15}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"$2"_"$1,\
"-hold_jid","S.09_FILTER_TO_SNV_ONLY_PASS_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".FILTER_TO_TITV_VCF_NOVEL.log",\
"'$SCRIPT_DIR'""/S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,$4,"'$DBSNP_129'""\n""sleep 3s"}'

### RUN TITV FOR THE PASSING SNVS THAT FALL IN UCSC CODING REGIONS THAT TOUCH EITHER THE BED OR TARGET FILE

## ALL SNVS TITV

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09-A.01-A.01_TITV_ALL_"$2"_"$1,\
"-hold_jid","S.09-A.01_FILTER_TO_SAMPLE_TITV_VCF_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".RUN_TITV_ALL.log",\
"'$SCRIPT_DIR'""/S.09-A.01-A.01_TITV_ALL.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

## ALL KNOWN SNVS TITV

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09-A.02-A.01_TITV_KNOWN_"$2"_"$1,\
"-hold_jid","S.09-A.02_FILTER_TO_SAMPLE_TITV_VCF_KNOWN_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".RUN_TITV_KNOWN.log",\
"'$SCRIPT_DIR'""/S.09-A.02-A.01_TITV_KNOWN.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

## ALL NOVEL SNVS TITV

awk 'BBEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.09-A.03-A.01_TITV_NOVEL_"$2"_"$1,\
"-hold_jid","S.09-A.03_FILTER_TO_SAMPLE_TITV_VCF_NOVEL_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".RUN_TITV_NOVEL.log",\
"'$SCRIPT_DIR'""/S.09-A.03-A.01_TITV_NOVEL.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

#####################################################################
################ CONVERT VCF FILES TO TABLES ########################
#####################################################################

## CONVERT INITIAL JOINT CALLED VCF TO TABLE##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","S.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_"$1,\
"-hold_jid","P.01_VARIANT_ANNOTATOR_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".VARIANT_TO_TABLE_COHORT_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.18_VARIANT_TO_TABLE_COHORT_ALL_SITES.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

## BGZIP INITIAL JOINT CALLED VCF TABLE##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","S.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES_"$1,\
"-hold_jid","S.18_VARIANT_TO_TABLE_COHORT_ALL_SITES_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1"\n""sleep 1s"}'

## TABIX INDEX INITIAL JOINT CALLED VCF TABLE##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 \
| uniq \
| awk '{print "qsub","-N","S.18-A.01-A.01_VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES_"$1,\
"-hold_jid","S.18-A.01_VARIANT_TO_TABLE_BGZIP_COHORT_ALL_SITES_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.18-A.01-A.01_VARIANT_TO_TABLE_TABIX_COHORT_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1"\n""sleep 1s"}'

######### SAMPLE ONLY ALL SITES FILE TO TABLE #################################

## CONVERT SAMPLE ONLY VCF TO TABLE##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.06-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_"$2"_"$1,\
"-hold_jid","S.06_FILTER_TO_SAMPLE_ALL_SITES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_SAMPLE_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.06-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES.sh",\
"'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

## BGZIP SAMPLE ONLY VCF TABLE##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.06-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES_"$2"_"$1,\
"-hold_jid","S.06-A.01_VARIANT_TO_TABLE_SAMPLE_ALL_SITES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.06-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

## TABIX INDEX SAMPLE ONLY VCF TABLE##

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| awk '{print "qsub","-N","S.06-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES_"$2"_"$1,\
"-hold_jid","S.06-A.01-A.01_VARIANT_TO_TABLE_BGZIP_SAMPLE_ALL_SITES_"$2"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES.log",\
"'$SCRIPT_DIR'""/S.06-A.01-A.01-A.01_VARIANT_TO_TABLE_TABIX_SAMPLE_ALL_SITES.sh",\
"'$TABIX_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

###################
##### ANNOVAR #####
###################

## RUN ANNOVAR

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{print "qsub","-N","S.07-A.01_RUN_ANNOVAR_"$2"_"$1,\
"-hold_jid","S.07_FILTER_TO_SAMPLE_VARIANTS_"$2"_"$1,\
"-pe slots 5",\
"-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".RUN_ANNOVAR.log",\
"'$SCRIPT_DIR'""/S.07-A.01_RUN_ANNOVAR.sh",\
"'$JAVA_1_6'","'$CIDRSEQSUITE_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

## REFORMAT ANNOVAR

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1,1 -k 2,2 \
| uniq \
| awk '{print "qsub","-N","S.07-A.01-A.01_REFORMAT_ANNOVAR_"$2"_"$1,\
"-hold_jid","S.07-A.01_RUN_ANNOVAR_"$2"_"$1,\
"-o","'$CORE_PATH'" "/" $1"/LOGS/"$2"_"$1".REFORMAT_ANNOVAR.log",\
"'$SCRIPT_DIR'""/S.07-A.01-A.01_REFORMAT_ANNOVAR.sh",\
"'$ANNOVAR_DIR'","'$CORE_PATH'",$1,$2"\n""sleep 3s"}'

##########################
##### QC REPORT PREP #####
##########################

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$20,$8,$21,$22,$23,$24}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 -k 3 \
| uniq \
| awk 'BEGIN {FS="\t"}
{print "qsub","-N","X.01-QC_REPORT_PREP_"$1"_"$3,\
"-hold_jid","S.07-A.01-A.01_REFORMAT_ANNOVAR_"$3"_"$1,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$3"_"$1".QC_REPORT_PREP.log",\
"'$SCRIPT_DIR'""/X.01-QC_REPORT_PREP.sh",\
"'$SAMTOOLS_DIR'","'$CORE_PATH'","'$DATAMASH_DIR'",$1,$2,$3,$4,$5,$6,$7"\n""sleep 30s"}'

### END PROJECT TASKS ###

awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
| sort -k 1 -k 2 \
| uniq \
| $DATAMASH_DIR/datamash -s -g 1 collapse 2 \
| awk 'BEGIN {FS="\t"}
gsub (/,/,",X.01-QC_REPORT_PREP_"$1"_",$2) \
{print "qsub","-N","X.01-X.01-END_PROJECT_TASKS_"$1,\
"-hold_jid","X.01-QC_REPORT_PREP_"$1"_"$2,\
"-o","'$CORE_PATH'/"$1"/LOGS/"$1".END_PROJECT_TASKS.log",\
"'$SCRIPT_DIR'""/X.01-X.01-END_PROJECT_TASKS.sh",\
"'$CORE_PATH'","'$DATAMASH_DIR'",$1"\n""sleep 3s"}'

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
