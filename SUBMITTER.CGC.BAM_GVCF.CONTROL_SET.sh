#!/usr/bin/env bash

# INPUT VARIABLES

	SAMPLE_SHEET=$1
	PED_FILE=$2
	QUEUE_LIST=$3 # optional. if no 3rd argument present then the default is cgc.q

		if [[ ! $QUEUE_LIST ]]
			then
			QUEUE_LIST="cgc.q"
		fi

	PRIORITY=$4 # optional. if no 4th argument present then the default is -15.
		# if you want to set this then you need to set the 3rd argument as well (even to the default)

			if [[ ! $PRIORITY ]]
				then
				PRIORITY="-15"
			fi

# CHANGE SCRIPT DIR TO WHERE YOU HAVE HAVE THE SCRIPTS BEING SUBMITTED

	SUBMITTER_SCRIPT_PATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

	SCRIPT_DIR="$SUBMITTER_SCRIPT_PATH/scripts"

##################
# CORE VARIABLES #
##################

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
	ALIGNMENT_CONTAINER="/mnt/clinical/ddl/NGS/CIDRSeqSuite/images/ddl_ce_control_align-0.0.1.simg"
	# contains the following software and is on Ubuntu 16.04.5 LTS
		# gatk 4.0.11.0 (base image). also contains the following.
			# Python 3.6.2 :: Continuum Analytics, Inc.
				# samtools 0.1.19
				# bcftools 0.1.19
				# bedtools v2.25.0
				# bgzip 1.2.1
				# tabix 1.2.1
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

##################
# PIPELINE FILES #
##################

	GENE_LIST="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/RefSeqGene.GRCh37.Ready.txt"
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
		$CORE_PATH/$PROJECT/REPORTS/{ALIGNMENT_SUMMARY,ANNOVAR,PICARD_DUPLICATES,TI_TV,VERIFYBAMID,VERIFYBAMID_CHR,ANEUPLOIDY_CHECK} \
		$CORE_PATH/$PROJECT/REPORTS/BAIT_BIAS/{METRICS,SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/PRE_ADAPTER/{METRICS,SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/BASECALL_Q_SCORE_DISTRIBUTION/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/REPORTS/BASE_DISTRIBUTION_BY_CYCLE/{METRICS,PDF} \
		$CORE_PATH/$PROJECT/REPORTS/COUNT_COVARIATES/GATK_REPORT \
		$CORE_PATH/$PROJECT/REPORTS/GC_BIAS/{METRICS,PDF,SUMMARY} \
		$CORE_PATH/$PROJECT/REPORTS/DEPTH_OF_COVERAGE/{TARGET,REFSEQ_CODING_PLUS_10bp} \
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
				$REF_GENOME
		}

	#######################################
	# run bqsr on the using bait bed file #
	#######################################

		RUN_BQSR ()
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

for SM_TAG in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
	do
		CREATE_SAMPLE_ARRAY
		FIX_BED_FILES
		echo sleep 0.1s
		RUN_BQSR
		echo sleep 0.1s
		APPLY_BQSR
		echo sleep 0.1s
		# SELECT_VERIFYBAMID_VCF
		# echo sleep 0.1s
		# RUN_VERIFYBAMID
		# echo sleep 0.1s
done

# ##### ALL H.00X SERIES OF SCRIPTS CAN BE RUN IN PARALLEL SINCE THEY ARE DEPENDENT ON FINAL BAM FILE GENERATION #####

# # Run Haplotype Caller in GVCF mode

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.01_HAPLOTYPE_CALLER_"$1"_"$2,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".HAPLOTYPE_CALLER.log",\
# "'$SCRIPT_DIR'""/H.01_HAPLOTYPE_CALLER.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # Run POST BQSR TABLE

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$19,$18}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
# print "qsub","-N","H.02_POST_BQSR_TABLE_"$2"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".POST_BQSR_TABLE.log",\
# "'$SCRIPT_DIR'""/H.02_POST_BQSR_TABLE.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3,INDEL[1],INDEL[2],$5"\n""sleep 1s"}'

# # Run ANALYZE COVARIATES

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($4,INDEL,";"); split($2,smtag,"[@-]"); \
# print "qsub","-N","H.02-A.01_ANALYZE_COVARIATES_"$2"_"$1,\
# "-hold_jid","H.02_POST_BQSR_TABLE_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".ANALYZE_COVARIATES.log",\
# "'$SCRIPT_DIR'""/H.02-A.01_ANALYZE_COVARIATES.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'",$1,$2,$3"\n""sleep 1s"}'

# # RUN DOC RefSeq CODING PLUS 10 BP FLANKS
# # This will specifically be for the RefSeg transcript IDs merged with UCSC canonical transcripts

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".DOC_CODING_10bpFLANKS.log",\
# "'$SCRIPT_DIR'""/H.03_DOC_CODING_10bpFLANKS.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3"\n""sleep 1s"}'

# # RUN ANEUPLOIDY_CHECK AFTER DOC RefSeq CODING PLUS 10 BP FLANKS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.03-A.01_DOC_CHROM_DEPTH_"$2"_"$1,\
# "-hold_jid","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".ANEUPLOIDY_CHECK.log",\
# "'$SCRIPT_DIR'""/H.03-A.01_CHROM_DEPTH.sh",\
# "'$CORE_PATH'","'$CYTOBAND_BED'","'$DATAMASH_DIR'","'$BEDTOOLS_DIR'",$1,$2"\n""sleep 1s"}'

# # RUN FORMATTING PER BASE COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{split($2,smtag,"[@]"); \
# print "qsub", "-N", "H.03-A.02_PER_BASE_" smtag[1] "_" smtag[2] "_" $1,\
# "-hold_jid", "H.03_DOC_CODING_10bpFLANKS_" $2 "_" $1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE.log",\
# "'$SCRIPT_DIR'""/H.03-A.02_PER_BASE.sh",\
# "'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2"\n""sleep 1s"}'

# # RUN FILTERING PER BASE COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{split($2,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02-A.01_PER_BASE_FILTER_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE_FILTER.log",\
# "'$SCRIPT_DIR'""/H.03-A.02-A.01_PER_BASE_FILTERED.sh",\
# "'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# # BGZIP PER BASE COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{split($2,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.02_PER_BASE_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE_BGZIP.log",\
# "'$SCRIPT_DIR'""/H.03-A.02-A.02_PER_BASE_BGZIP.sh",\
# "'$CORE_PATH'","'$TABIX_DIR'",$1,$2"\n""sleep 1s"}'

# # TABIX PER BASE COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{split($2,smtag,"[@]"); \
# print "qsub","-N","H.03-A.02-A.02-A.01_PER_BASE_TABIX_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.02-A.02_PER_BASE_BGZIP_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_BASE_TABIX.log",\
# "'$SCRIPT_DIR'""/H.03-A.02-A.02-A.01_PER_BASE_TABIX.sh",\
# "'$CORE_PATH'","'$TABIX_DIR'",$1,$2"\n""sleep 1s"}'

# # RUN FORMATTING PER CODING INTERVAL COVERAGE WITH GENE NAME ANNNOTATION

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{split($2,smtag,"[@]"); \
# print "qsub","-N","H.03-A.03_PER_INTERVAL_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03_DOC_CODING_10bpFLANKS_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_INTERVAL.log",\
# "'$SCRIPT_DIR'""/H.03-A.03_PER_INTERVAL.sh",\
# "'$CORE_PATH'","'$BEDTOOLS_DIR'","'$CODING_BED'",$1,$2"\n""sleep 1s"}'

# # RUN FILTERING PER CODING INTERVAL COVERAGE WITH GENE NAME ANNNOTATION WITH LESS THAN 30x

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1,1 -k 2,2 \
# | uniq \
# | awk '{split($2,smtag,"[@]"); \
# print "qsub","-N","H.03-A.03_PER_INTERVAL_FILTER_"smtag[1]"_"smtag[2]"_"$1,\
# "-hold_jid","H.03-A.03_PER_INTERVAL_"smtag[1]"_"smtag[2]"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".PER_INTERVAL_FILTER.log",\
# "'$SCRIPT_DIR'""/H.03-A.03-A.01_PER_INTERVAL_FILTERED.sh",\
# "'$CORE_PATH'",$1,$2"\n""sleep 1s"}'

# # RUN DOC TARGET BED (Generally this with all RefGene coding exons unless it becomes targeted)

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.05_DOC_TARGET_BED_"$2"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".DOC_TARGET_BED.log",\
# "'$SCRIPT_DIR'""/H.05_DOC_TARGET_BED.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$GENE_LIST'",$1,$2,$3"\n""sleep 1s"}'

# # RUN COLLECT MULTIPLE METRICS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$18,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.06_COLLECT_MULTIPLE_METRICS_"$2"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".COLLECT_MULTIPLE_METRICS.log",\
# "'$SCRIPT_DIR'""/H.06_COLLECT_MULTIPLE_METRICS.sh",\
# "'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3,$4,$5"\n""sleep 1s"}'

# # RUN COLLECT HS METRICS

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.07_COLLECT_HS_METRICS_"$2"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".COLLECT_HS_METRICS.log",\
# "'$SCRIPT_DIR'""/H.07_COLLECT_HS_METRICS.sh",\
# "'$JAVA_1_8'","'$PICARD_DIR'","'$CORE_PATH'","'$SAMTOOLS_DIR'",$1,$2,$3"\n""sleep 1s"}'

# # RUN SELECT VERIFYBAM ID VCF

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8,$12,$15}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.08_SELECT_VERIFYBAMID_VCF_"$2"_"$1,\
# "-hold_jid","G.01_FINAL_BAM_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".SELECT_VERIFYBAMID_VCF.log",\
# "'$SCRIPT_DIR'""/H.08_SELECT_VERIFYBAMID_VCF.sh",\
# "'$JAVA_1_8'","'$GATK_DIR'","'$CORE_PATH'","'$VERIFY_VCF'",$1,$2,$3,$4"\n""sleep 1s"}'

# # RUN VERIFYBAMID ALL

# awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,$8}' \
# ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt \
# | sort -k 1 -k 2 \
# | uniq \
# | awk '{split($2,smtag,"[@-]"); \
# print "qsub","-N","H.08-A.01_VERIFYBAMID_"$2"_"$1,\
# "-hold_jid","H.08_SELECT_VERIFYBAMID_VCF_"$2"_"$1,\
# "-o","'$CORE_PATH'/"$1"/LOGS/"$2"_"$1".VERIFYBAMID.log",\
# "'$SCRIPT_DIR'""/H.08-A.01_VERIFYBAMID.sh",\
# "'$CORE_PATH'","'$VERIFY_DIR'",$1,$2"\n""sleep 1s"}'

# ###################################################
# ### RUN VERIFYBAM ID PER CHROMOSOME - VITO ########
# ###################################################

# CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM ()
# {
# SAMPLE_INFO_ARRAY_VERIFY_BAM=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$20,$8,$12,$15}' ~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)
# }

# CALL_SELECT_VERIFY_BAM ()
# {
# echo \
# qsub \
# -N H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
# -hold_jid G.01_FINAL_BAM_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.SELECT_VERIFYBAMID_chr$CHROMOSOME.log \
# $SCRIPT_DIR/H.09_SELECT_VERIFYBAMID_VCF_CHR.sh \
# $JAVA_1_8 $GATK_DIR $CORE_PATH $VERIFY_VCF \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[3]} \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[4]} $CHROMOSOME
# }

# CALL_VERIFYBAMID ()
# {
# echo \
# qsub \
# -N H.09-A.01_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
# -hold_jid H.09_SELECT_VERIFYBAMID_VCF_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}_chr$CHROMOSOME \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.VERIFYBAMID_chr$CHROMOSOME.log \
# $SCRIPT_DIR/H.09-A.01_VERIFYBAMID_CHR.sh \
# $CORE_PATH $VERIFY_DIR \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]} \
# $CHROMOSOME
# }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
# do
# CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
# 	for CHROMOSOME in {1..22}
# 		do
# 		CALL_SELECT_VERIFY_BAM
# 		echo sleep 1s
# 		CALL_VERIFYBAMID
# 		echo sleep 1s
# 	done
# done

# #####################################################
# ### JOIN THE PER CHROMOSOME VERIFYBAMID REPORTS #####
# #####################################################

# BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR ()
# {
# 	for PROJECT in $(awk 'BEGIN {FS=","} NR>1 {print $1}' $SAMPLE_SHEET | sort | uniq )
# 	do
# 	HOLD_ID_PATH="-hold_jid "
# 	for CHROMOSOME in {{1..22},{X,Y}};
#  	do
#  		HOLD_ID_PATH=$HOLD_ID_PATH"H.09-A.01_VERIFYBAMID_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}"_"${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}"_"chr$CHROMOSOME","
#  	done
#  done
# }

#  CAT_VERIFYBAMID_CHR ()
#  {
# echo \
# qsub \
# -N H.09-A.01-A.01_JOIN_VERIFYBAMID_${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} \
# $HOLD_ID_PATH \
# -o $CORE_PATH/${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}/LOGS/${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}_${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]}.CAT_VERIFYBAMID_CHR.log \
# $SCRIPT_DIR/H.09-A.01-A.01_CAT_VERIFYBAMID_CHR.sh \
# $CORE_PATH \
# ${SAMPLE_INFO_ARRAY_VERIFY_BAM[0]} ${SAMPLE_INFO_ARRAY_VERIFY_BAM[2]}
#  }

# for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
#  do
#  	CREATE_SAMPLE_INFO_ARRAY_VERIFY_BAM
# 	BUILD_HOLD_ID_PATH_CAT_VERIFYBAMID_CHR
# 	CAT_VERIFYBAMID_CHR
# 	echo sleep 1s
#  done

#############################################
