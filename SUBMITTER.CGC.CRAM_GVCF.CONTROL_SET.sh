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

	SCRIPT_DIR="$SUBMITTER_SCRIPT_PATH/scripts_cram_gvcf"

##################
# CORE VARIABLES #
##################

	# GVCF PAD. CURRENTLY KEEPING THIS AS A STATIC VARIABLE

		GVCF_PAD="250"

	## This will always put the current working directory in front of any directory for PATH
	## added /bin for RHEL6

		export PATH=".:$PATH:/bin"

	# where the input/output sequencing data will be located.
		
		CORE_PATH="/mnt/clinical/ddl/NGS/Exome_Data"

	# Directory where NovaSeqa runs are located.

		NOVASEQ_REPO="/mnt/instrument_files/novaseq"

	# this is just for note taking where the control data set will reside. not used in this pipeline.

		CONTROL_REPO="/mnt/clinical/ddl/NGS/Exome_Resources/CONTROL_REPO_TWIST"

	# used for tracking in the read group header of the cram file

		PIPELINE_VERSION=`git --git-dir=$SCRIPT_DIR/../.git --work-tree=$SCRIPT_DIR/.. log --pretty=format:'%h' -n 1`

	# load gcc for programs like verifyBamID
	## this will get pushed out to all of the compute nodes since I specify env var to pushed out with qsub
		
		module load gcc/7.2.0

	# explicitly setting this b/c not everybody has had the $HOME directory transferred and I'm not going to through
	# and figure out who does and does not have this set correctly
		
		umask 0007

	# SUBMIT TIMESTAMP

		SUBMIT_STAMP=`date '+%s'`

	# SUBMITTER_ID
		
		SUBMITTER_ID=`whoami`

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
		| join -1 8 -2 2 -e '-'  -t $'\t' \
		-o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.1,2.3,2.4,2.5,2.6' \
		~/CGC_PIPELINE_TEMP/SORTED.$MANIFEST_PREFIX.txt /dev/stdin \
		>| ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt
	}

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
#################################

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

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
		do
			SETUP_PROJECT
done

################################
##### CRAM FILE GENERATION #####
###############################################################################################
##### NOTE: THE CRAM FILE IS THE END PRODUCT BUT THE BAM FILE IS USED FOR OTHER PROCESSES #####
##### SOME PROGRAMS CAN'T TAKE IN CRAM AS AN INPUT ############################################
###############################################################################################

	########################################################################################
	# create an array at the platform level so that bwa mem can add metadata to the header #
	########################################################################################

		CREATE_PLATFORM_UNIT_ARRAY ()
		{
			PLATFORM_UNIT_ARRAY=(`awk 1 ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
			| awk 'BEGIN {FS="\t"} $8$2$3$4=="'$PLATFORM_UNIT'" {split($19,INDEL,";"); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2],$$20,$21,$22,$23,$24}' \
			| sort \
			| uniq`)

				#  1  Project=the Seq Proj folder name
				
					PROJECT=${PLATFORM_UNIT_ARRAY[0]}

				#  2  FCID=flowcell that sample read group was performed on
				
					FCID=${PLATFORM_UNIT_ARRAY[1]}

				#  3  Lane=lane of flowcell that sample read group was performed on

					LANE=${PLATFORM_UNIT_ARRAY[2]}

				#  4  Index=sample barcode

					INDEX=${PLATFORM_UNIT_ARRAY[3]}

				#  5  Platform=type of sequencing chemistry matching SAM specification

					PLATFORM=${PLATFORM_UNIT_ARRAY[4]}

				#  6  Library_Name=library group of the sample read group, 
					# Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not

					LIBRARY=${PLATFORM_UNIT_ARRAY[5]}

				#  7  Date=should be the run set up date, but doesn't have to be

					RUN_DATE=${PLATFORM_UNIT_ARRAY[6]}

				#  8  SM_Tag=sample ID

					SM_TAG=${PLATFORM_UNIT_ARRAY[7]}
						
						# sge sm tag. If there is an @ in the qsub or holdId name it breaks

							SGE_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g')

				#  9  Center=the center/funding mechanism

					CENTER=${PLATFORM_UNIT_ARRAY[8]}

				# 10  Description=Generally we use to denote the sequencer setting (e.g. rapid run)
				# “HiSeq-X”, “HiSeq-4000”, “HiSeq-2500”, “HiSeq-2000”, “NextSeq-500”, or “MiSeq”.

					SEQUENCER_MODEL=${PLATFORM_UNIT_ARRAY[9]}

				########################
				# 11  Seq_Exp_ID: SKIP #
				########################

				# 12  Genome_Ref=the reference genome used in the analysis pipeline

					REF_GENOME=${PLATFORM_UNIT_ARRAY[10]}

				#####################################
				# 13  Operator: SKIP ################
				# 14  Extra_VCF_Filter_Params: SKIP #
				#####################################

				# 15  TS_TV_BED_File=refseq (select) cds plus other odds and ends (.e.g. missing omim))

					TITV_BED=${PLATFORM_UNIT_ARRAY[11]}

				# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
				# Used for limited where to run base quality score recalibration on where to create gvcf files.
				
					BAIT_BED=${PLATFORM_UNIT_ARRAY[12]}

				# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.
				
					TARGET_BED=${PLATFORM_UNIT_ARRAY[13]}

				# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in base call quality score recalibration.

					DBSNP=${PLATFORM_UNIT_ARRAY[14]}

				# 19  KNOWN_INDEL_FILES=used for BQSR masking

					KNOWN_INDEL_1=${PLATFORM_UNIT_ARRAY[15]}
					KNOWN_INDEL_2=${PLATFORM_UNIT_ARRAY[16]}

				# 20 FAMILY

					FAMILY=${PLATFORM_UNIT_ARRAY[17]}

				# 21 MOM

					MOM=${PLATFORM_UNIT_ARRAY[17]}

				# 22 DAD

					DAD=${PLATFORM_UNIT_ARRAY[17]}

				# 23 GENDER

					GENDER=${PLATFORM_UNIT_ARRAY[17]}

				# 24 PHENOTYPE

					PHENOTYPE=${PLATFORM_UNIT_ARRAY[17]}
		}

	########################################################################
	### Use bwa mem to do the alignments; ##################################
	### pipe to samblaster to add mate tags; ###############################
	### pipe to picard's AddOrReplaceReadGroups to handle the bam header ###
	########################################################################

		RUN_BWA ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N A.01-BWA"_"$SGE_SM_TAG"_"$FCID"_"$LANE"_"$INDEX \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"_"$FCID"_"$LANE"_"$INDEX"-BWA.log" \
			$SCRIPT_DIR/A.01_BWA.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$FCID \
				$LANE \
				$INDEX \
				$PLATFORM \
				$LIBRARY \
				$RUN_DATE \
				$SM_TAG \
				$CENTER \
				$SEQUENCER_MODEL \
				$REF_GENOME \
				$PIPELINE_VERSION \
				$BAIT_BED \
				$TARGET_BED \
				$TITV_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP \
				$NOVASEQ_REPO
		}

		for PLATFORM_UNIT in $(awk 'BEGIN {FS=","} NR>1 {print $8$2$3$4}' $SAMPLE_SHEET | sort | uniq );
			do
				CREATE_PLATFORM_UNIT_ARRAY
				mkdir -p $CORE_PATH/$PROJECT/LOGS/$SM_TAG
				RUN_BWA
				echo sleep 0.1s
		done

	#########################################################################################
	# Merge files and mark duplicates using picard duplictes with queryname sorting #########
	# do coordinate sorting with sambamba ###################################################
	#########################################################################################
	# I am setting the heap space and garbage collector threads for picard now now ##########
	# doing this does drastically decrease the load average ( the gc thread specification ) #
	#########################################################################################
	# create a hold job id qsub command line based on the number of #########################
	# submit merging the bam files created by bwa mem above #################################
	# only launch when every lane for a sample is done being processed by bwa mem ###########
	# I want to clean this up eventually and get away from using awk to print the qsub line #
	#########################################################################################

		awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'BEGIN {FS=","; OFS="\t"} NR>1 {print $1,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam",$8,$10}' \
			| awk 'BEGIN {OFS="\t"} {sub(/@/,"_",$5)} {print $1,$2,$3,$4,$5,$6}' \
			| sort -k 1,1 -k 2,2 -k 3,3 -k 6,6 \
			| uniq \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				-s \
				-g 1,2 \
				collapse 3 \
				collapse 4 \
				unique 5 \
				unique 6 \
			| awk 'BEGIN {FS="\t"} \
				gsub(/,/,",A.01-BWA_"$5"_",$3) \
				gsub(/,/,",INPUT=" "'$CORE_PATH'" "/" $1"/TEMP/",$4) \
				{print "qsub",\
				"-S /bin/bash",\
				"-cwd",\
				"-V",\
				"-v SINGULARITY_BINDPATH=/mnt:/mnt",\
				"-q","'$QUEUE_LIST'",\
				"-p","'$PRIORITY'",\
				"-j y",\
				"-N","B.01-MARK_DUPLICATES_"$5"_"$1,\
				"-o","'$CORE_PATH'/"$1"/LOGS/"$2"/"$2"-MARK_DUPLICATES.log",\
				"-hold_jid","A.01-BWA_"$5"_"$3, \
				"'$SCRIPT_DIR'""/B.01_MARK_DUPLICATES.sh",\
				"'$ALIGNMENT_CONTAINER'",\
				"'$CORE_PATH'",\
				$1,\
				$2,\
				"'$SAMPLE_SHEET'",\
				"'$SUBMIT_STAMP'",\
				$6,\
				"INPUT=" "'$CORE_PATH'" "/" $1"/TEMP/"$4"\n""sleep 0.1s"}'

	###############################################
	# fix common formatting problems in bed files #
	# merge bait to target for gvcf creation, pad #
	# create picard style interval files ##########
	###############################################

		FIX_BED_FILES ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-FIX_BED_FILES.log" \
			-hold_jid B.01-MARK_DUPLICATES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/C.01_FIX_BED.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$TARGET_BED \
				$BAIT_BED \
				$TITV_BED \
				$CYTOBAND_BED \
				$REF_GENOME \
				$PADDING_LENGTH \
				$GVCF_PAD
		}

	#######################################
	# run bqsr on the using bait bed file #
	#######################################

		PERFORM_BQSR ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N D.01-PERFORM_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-PERFORM_BQSR.log" \
			-hold_jid B.01-MARK_DUPLICATES"_"$SGE_SM_TAG"_"$PROJECT,C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/D.01_PERFORM_BQSR.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$KNOWN_INDEL_1 \
				$KNOWN_INDEL_2 \
				$DBSNP \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################
	# use a 4 bin q score scheme #
	# remove indel Q scores ######
	# retain original Q score  ###
	##############################

		APPLY_BQSR ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-APPLY_BQSR.log" \
			-hold_jid D.01-PERFORM_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/E.01_APPLY_BQSR.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#####################################################
	# create a lossless cram, although the bam is lossy #
	#####################################################

		BAM_TO_CRAM ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-BAM_TO_CRAM.log" \
			-hold_jid E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/F.01_BAM_TO_CRAM.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	# ##########################################################################################
	# # index the cram file and copy it so that there are both *crai and cram.crai *extensions #
	# ##########################################################################################

		INDEX_CRAM ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N G.01-INDEX_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-INDEX_CRAM.log" \
			-hold_jid F.01-BAM_TO_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/G.01_INDEX_CRAM.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

for SM_TAG in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
	do
		CREATE_SAMPLE_ARRAY
		FIX_BED_FILES
		echo sleep 0.1s
		PERFORM_BQSR
		echo sleep 0.1s
		APPLY_BQSR
		echo sleep 0.1s
		BAM_TO_CRAM
		echo sleep 0.1s
		INDEX_CRAM
done

########################################################################################
##### BAM/CRAM FILE RELATED METRICS ####################################################
##### NOTE: SOME PROGRAMS CAN ONLY BE RAN ON THE BAM FILE AND NOT ON THE CRAM FILE #####
##### I WILL COMMENT ON WHICH IS WHICH #################################################
########################################################################################

	################################################################################
	# COLLECT MULTIPLE METRICS  ####################################################
	# again used bait bed file here instead of target b/c target could be anything #
	# ti/tv bed is unrelated to the capture really #################################
	# uses the CRAM file as the input ##############################################
	################################################################################

		COLLECT_MULTIPLE_METRICS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.01-COLLECT_MULTIPLE_METRICS"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-COLLECT_MULTIPLE_METRICS.log" \
			-hold_jid G.01-INDEX_CRAM"_"$SGE_SM_TAG"_"$PROJECT,C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.01_COLLECT_MULTIPLE_METRICS.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$DBSNP \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#########################################
	# COLLECT HS METRICS  ###################
	# bait bed is the bait bed file #########
	# titv bed files is the target bed file #
	# uses the CRAM file as the input #######
	#########################################

		COLLECT_HS_METRICS ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.02-COLLECT_HS_METRICS"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-COLLECT_HS_METRICS.log" \
			-hold_jid G.01-INDEX_CRAM"_"$SGE_SM_TAG"_"$PROJECT,C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.02_COLLECT_HS_METRICS.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$BAIT_BED \
				$TITV_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################################################################
	# CREATE DEPTH OF COVERAGE FOR TARGET BED PADDED WITH THE INPUT FROM THE GUI #
	# uses a gatk 3.7 container ##################################################
	# input is the BAM file #################################################################################
	# Generally this with all RefSeq Select CDS exons + missing OMIM unless it becomes targeted, e.g a zoom #
	# uses the BAM file as the input ########################################################################
	#########################################################################################################

		DOC_TARGET ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.03-DOC_TARGET"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-DOC_TARGET.log" \
			-hold_jid G.01-INDEX_CRAM"_"$SGE_SM_TAG"_"$PROJECT,C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.03_DOC_TARGET_PADDED_BED.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$TARGET_BED \
				$PADDING_LENGTH \
				$GENE_LIST \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#################################################################################################
	# CREATE VCF FOR VERIFYBAMID METRICS ############################################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	#################################################################################################

		SELECT_VERIFYBAMID_VCF ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.04-SELECT_VERIFYBAMID_VCF"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-SELECT_VERIFYBAMID_VCF.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.04_SELECT_VERIFYBAMID_VCF.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$VERIFY_VCF \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	###############################
	# RUN VERIFYBAMID #############
	# THIS RUNS OFF OF A BAM FILE #
	###############################

		RUN_VERIFYBAMID ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.04-A.01-RUN_VERIFYBAMID"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-VERIFYBAMID.log" \
			-hold_jid H.04-SELECT_VERIFYBAMID_VCF"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.04-A.01_VERIFYBAMID.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##############################################################################
	# CREATE DEPTH OF COVERAGE FOR CODING BED PADDED WITH THE INPUT FROM THE GUI #
	# uses a gatk 3.7 container ##################################################
	# input is the BAM file ######################################################
	# This with all RefSeq Select CDS exons + missing OMIM, etc. #################
	# uses the BAM file as the input #############################################
	##############################################################################

		DOC_CODING ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-DOC_CODING.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,G.01-INDEX_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05_DOC_CODING_PADDED.sh \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$CODING_BED \
				$PADDING_LENGTH \
				$GENE_LIST \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#########################################################
	# DO AN ANEUPLOIDY CHECK ON TARGET BED FILE DOC OUTPUT  #
	#########################################################

		ANEUPLOIDY_CHECK ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.01_CHROM_DEPTH"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-ANEUPLOIDY_CHECK.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.01_CHROM_DEPTH.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	########################################################################################
	# FORMATTING PER BASE COVERAGE AND ADDING GENE NAME, TRANSCRIPT, EXON, ETC ANNNOTATION #
	########################################################################################

		ANNOTATE_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02_ANNOTATE_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-ANNOTATE_PER_BASE.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02_ANNOTATE_PER_BASE.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	##########################################################################
	# FILTER PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x #
	##########################################################################

		FILTER_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02-A.01_FILTER_ANNOTATED_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-FILTER_ANNOTATED_PER_BASE.log" \
			-hold_jid H.05-A.02_ANNOTATE_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02-A.01_FILTER_ANNOTATED_PER_BASE.sh \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	######################################################
	# BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION #
	######################################################

		BGZIP_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02-A.02_BGZIP_ANNOTATED_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-BGZIP_ANNOTATED_PER_BASE.log" \
			-hold_jid H.05-A.02_ANNOTATE_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02-A.02_BGZIP_ANNOTATED_PER_BASE.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	######################################################
	# TABIX PER BASE COVERAGE WITH GENE NAME ANNNOTATION #
	######################################################

		TABIX_ANNOTATED_PER_BASE_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.02-A.02-A.01_TABIX_ANNOTATED_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-TABIX_ANNOTATED_PER_BASE.log" \
			-hold_jid H.05-A.02-A.02_BGZIP_ANNOTATED_PER_BASE"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.02-A.02-A.01_TABIX_ANNOTATED_PER_BASE.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	###################################################################################################
	# FORMATTING PER CODING INTERVAL COVERAGE AND ADDING GENE NAME, TRANSCRIPT, EXON, ETC ANNNOTATION #
	###################################################################################################

		ANNOTATE_PER_INTERVAL_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.03_ANNOTATE_PER_INTERVAL"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-ANNOTATE_PER_INTERVAL.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,H.05-DOC_CODING"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.03_ANNOTATE_PER_INTERVAL.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	##################################################################################################
	# FILTER ANNOTATED PER CODING INTERVAL COVERAGE TO INTERVALS WHERE LESS 100% OF BASES ARE AT 30X #
	##################################################################################################

		FILTER_ANNOTATED_PER_INTERVAL_REPORT ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.05-A.03-A.01_FILTER_ANNOTATED_PER_INTERVAL"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-FILTER_ANNOTATED_PER_INTERVAL.log" \
			-hold_jid H.05-A.03_ANNOTATE_PER_INTERVAL"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.05-A.03-A.01_FILTER_ANNOTATED_PER_INTERVAL.sh \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$CODING_BED \
				$PADDING_LENGTH
		}

	#################################################################################################
	# CREATE VCF PER CHROMOSOME AND RUN VERIFYBAMID ON THEM ######################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	# SCRIPT READS BAIT BED FILE, GRABS THE CHROMOSOMES AND RUNS A FOR LOOP FOR BOTH THINGS #########
	# USES BAM FILE AS THE INPUT ####################################################################
	#################################################################################################

		VERIFYBAMID_PER_AUTOSOME ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.06-SELECT_VERIFYBAMID_PER_AUTOSOME"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-SELECT_VERIFYBAMID_PER_AUTOSOME.log" \
			-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.06_VERIFYBAMID_PER_AUTO.sh \
				$ALIGNMENT_CONTAINER \
				$GATK_3_7_0_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$VERIFY_VCF \
				$BAIT_BED \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	#################################################################################################
	# CREATE VCF PER CHROMOSOME AND RUN VERIFYBAMID ON THEM ######################################
	# USE THE BAIT BED FILE #########################################################################
	# THE TARGET BED COULD BE MODIFIED TO BE TOO SMALL TO BE USEFUL HERE ############################
	# TI/TV BED FILE HAS TOO MUCH UNCERTAINTY SINCE IT DOES NOT HAE ANYTHING TO DO WITH THE CAPTURE #
	# SCRIPT READS BAIT BED FILE, GRABS THE CHROMOSOMES AND RUNS A FOR LOOP FOR BOTH THINGS #########
	# USES BAM FILE AS THE INPUT ####################################################################
	#################################################################################################

		CAT_VERIFYBAMID_PER_AUTOSOME ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.06-A.01-CAT_VERIFYBAMID_AUTOSOME"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-CAT_VERIFYBAMID_AUTOSOME.log" \
			-hold_jid H.06-SELECT_VERIFYBAMID_PER_AUTOSOME"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.06-A.01_CAT_VERIFYBAMID_AUTO.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$BAIT_BED
		}

for SM_TAG in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
	do
		CREATE_SAMPLE_ARRAY
		COLLECT_MULTIPLE_METRICS
		echo sleep 0.1s
		COLLECT_HS_METRICS
		echo sleep 0.1s
		DOC_TARGET
		echo sleep 0.1s
		SELECT_VERIFYBAMID_VCF
		echo sleep 0.1s
		RUN_VERIFYBAMID
		echo sleep 0.1s
		DOC_CODING
		echo sleep 0.1s
		ANEUPLOIDY_CHECK
		echo sleep 0.1s
		ANNOTATE_PER_BASE_REPORT
		echo sleep 0.1s
		FILTER_ANNOTATED_PER_BASE_REPORT
		echo sleep 0.1s
		BGZIP_ANNOTATED_PER_BASE_REPORT
		echo sleep 0.1s
		TABIX_ANNOTATED_PER_BASE_REPORT
		echo sleep 0.1s
		ANNOTATE_PER_INTERVAL_REPORT
		echo sleep 0.1s
		FILTER_ANNOTATED_PER_INTERVAL_REPORT
		echo sleep 0.1s
		VERIFYBAMID_PER_AUTOSOME
		echo sleep 0.1s
		CAT_VERIFYBAMID_PER_AUTOSOME
		echo sleep 0.1s
done

#####################################################################
# HAPLOTYPE CALLER SCATTER ##########################################
#####################################################################
# INPUT IS THE BAM FILE #############################################
# THE BED FILE FOR THE GVCF INTERVALS IS ############################
# THE BAIT BED PLUS THE CODING BED FILE CONCATENTATED TOGETHER ######
# THEN A 250 BP PAD ADDED AND THEN MERGED FOR OVERLAPPING INTERVALS #
#####################################################################

	CALL_HAPLOTYPE_CALLER ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N H.07-HAPLOTYPE_CALLER"_"$SGE_SM_TAG"_"$PROJECT"_chr"$CHROMOSOME \
			-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-HAPLOTYPE_CALLER_chr"$CHROMOSOME".log" \
		-hold_jid C.01-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT,E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
		$SCRIPT_DIR/H.07_HAPLOTYPE_CALLER_SCATTER.sh \
			$GATK_3_7_0_CONTAINER \
			$CORE_PATH \
			$PROJECT \
			$SM_TAG \
			$REF_GENOME \
			$CODING_BED \
			$BAIT_BED \
			$CHROMOSOME \
			$GVCF_PAD \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

# Take the samples bait bed file, create a list of unique chromosome to use as a scatter for haplotype_caller_scatter

for SM_TAG in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
	do
	CREATE_SAMPLE_ARRAY
		for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
			| sed -r 's/[[:space:]]+/\t/g' \
			| sed 's/chr//g' \
			| grep -v "MT" \
			| cut -f 1 \
			| sort \
			| uniq \
			| singularity exec $ALIGNMENT_CONTAINER datamash \
				collapse 1 \
			| sed 's/,/ /g');
			do
				CALL_HAPLOTYPE_CALLER
				echo sleep 0.1s
		done
done

###########################
# HAPLOTYPE CALLER GATHER #
################################################################################
# GATHER UP THE PER SAMPLE PER CHROMOSOME GVCF FILES INTO A SINGLE SAMPLE GVCF #
################################################################################

	BUILD_HOLD_ID_PATH ()
	{
		HOLD_ID_PATH="-hold_jid "

		for CHROMOSOME in $(sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' $BAIT_BED \
								| sed -r 's/[[:space:]]+/\t/g' \
								| cut -f 1 \
								| sed 's/chr//g' \
								| grep -v "MT" \
								| sort \
								| uniq \
								| singularity exec $ALIGNMENT_CONTAINER datamash \
									collapse 1 \
								| sed 's/,/ /g');
			do
				HOLD_ID_PATH=$HOLD_ID_PATH"H.07-HAPLOTYPE_CALLER_"$SM_TAG"_"$PROJECT"_chr"$CHROMOSOME","
				HOLD_ID_PATH=`echo $HOLD_ID_PATH | sed 's/@/_/g'`
		done
	}

	CALL_HAPLOTYPE_CALLER_GVCF_GATHER ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N H.01-A.01_HAPLOTYPE_CALLER_GVCF_GATHER"_"$SGE_SM_TAG"_"$PROJECT \
			-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG-HAPLOTYPE_CALLER_GVCF_GATHER.log \
		${HOLD_ID_PATH} \
		$SCRIPT_DIR/H.07-A.01_HAPLOTYPE_CALLER_GVCF_GATHER.sh \
			$GATK_3_7_0_CONTAINER \
			$CORE_PATH \
			$PROJECT \
			$SM_TAG \
			$REF_GENOME \
			$BAIT_BED \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

	CALL_HAPLOTYPE_CALLER_BAM_GATHER ()
	{
		echo \
		qsub \
			$QSUB_ARGS \
		-N H.07-A.02_HAPLOTYPE_CALLER_BAM_GATHER"_"$SGE_SM_TAG"_"$PROJECT \
			-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG-HAPLOTYPE_CALLER_BAM_GATHER.log \
		${HOLD_ID_PATH} \
		$SCRIPT_DIR/H.07-A.02_HAPLOTYPE_CALLER_BAM_GATHER.sh \
			$ALIGNMENT_CONTAINER \
			$CORE_PATH \
			$PROJECT \
			$SM_TAG \
			$BAIT_BED \
			$SAMPLE_SHEET \
			$SUBMIT_STAMP
	}

	########################################################
	# create a lossless HC cram, although the bam is lossy #
	########################################################

		HC_BAM_TO_CRAM ()
		{
			echo \
			qsub \
				$QSUB_ARGS \
			-N H.07-A.02-A.01_HAPLOTYPE_CALLER_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-HC_BAM_TO_CRAM.log" \
			-hold_jid H.07-A.02_HAPLOTYPE_CALLER_BAM_GATHER"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.07-A.02-A.01_HAPLOTYPE_CALLER_CRAM.sh \
				$ALIGNMENT_CONTAINER \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

	##########################################################################################
	# index the cram file and copy it so that there are both *crai and cram.crai *extensions #
	##########################################################################################

		HC_INDEX_CRAM ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
			-N H.01-A.02-A.01-A.01_INDEX_HAPLOTYPE_CALLER_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-HC_INDEX_CRAM.log" \
				-j y \
			-hold_jid H.01-A.02-A.01_HAPLOTYPE_CALLER_CRAM"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/H.01-A.02-A.01-A.01_INDEX_HAPLOTYPE_CALLER_CRAM.sh \
				$SAMTOOLS_DIR \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME
		}

for SM_TAG in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
	do
		CREATE_SAMPLE_ARRAY
		BUILD_HOLD_ID_PATH
		CALL_HAPLOTYPE_CALLER_GVCF_GATHER
		echo sleep 0.1s
		CALL_HAPLOTYPE_CALLER_BAM_GATHER
		echo sleep 0.1s
		HC_BAM_TO_CRAM
		echo sleep 0.1s
done
