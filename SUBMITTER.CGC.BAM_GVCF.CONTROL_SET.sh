#!/usr/bin/env bash

# INPUT VARIABLES

	SAMPLE_SHEET=$1
	PED_FILE=$2
	PRIORITY=$3 # optional. if no 2nd argument present then the default is -15

		# if there is no 3rd argument present then use the number for priority
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

	# where the input/output sequencing data will be located.
		
		CORE_PATH="/mnt/clinical/ddl/NGS/Exome_Data"

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

	# the sqe queue(s) that you want to submit to

		QUEUE_LIST="cgc.q"

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
	CYTOBAND_BED="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/GRCh37.Cytobands.bed"
	HAPMAP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/hapmap_3.3.b37.vcf"
	OMNI_1KG="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_omni2.5.b37.vcf"
	HI_CONF_1KG_PHASE1_SNP="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/1000G_phase1.snps.high_confidence.b37.vcf"
	MILLS_1KG_GOLD_INDEL="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/Mills_and_1000G_gold_standard.indels.b37.vcf"
	PHASE3_1KG_AUTOSOMES="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
	DBSNP_129="/mnt/clinical/ddl/NGS/Exome_Resources/PIPELINE_FILES/dbsnp_138.b37.excluding_sites_after_129.vcf"

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

	CREATE_SAMPLE_INFO_ARRAY ()
	{
		SAMPLE_INFO_ARRAY=(`awk 'BEGIN {FS="\t"; OFS="\t"} $8=="'$SAMPLE'" {print $1,$20,$8,$15,$16}' \
			~/CGC_PIPELINE_TEMP/$MANIFEST_PREFIX.$PED_PREFIX.join.txt`)

		PROJECT=${SAMPLE_INFO_ARRAY[0]}
		FAMILY=${SAMPLE_INFO_ARRAY[1]}
		SM_TAG=${SAMPLE_INFO_ARRAY[2]}
		BAIT_BED_FILE=${SAMPLE_INFO_ARRAY[3]}
		TARGET_BED_FILE=${SAMPLE_INFO_ARRAY[4]}
	}

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

#####################################################################################################
# PAD THE REFSEQ canonical transcript bed file by 10 bases. #########################################
# can make this as an input variable with a default value 10 if i have to ever give more than 0 effs.
############# this should be it's own job ###########################################################
#####################################################################################################

	SETUP_PROJECT ()
	{
		FORMAT_MANIFEST
		MERGE_PED_MANIFEST
		CREATE_SAMPLE_INFO_ARRAY
		MAKE_PROJ_DIR_TREE
	}

for SAMPLE in $(awk 'BEGIN {FS=","} NR>1 {print $8}' $SAMPLE_SHEET | sort | uniq );
do
SETUP_PROJECT
done

############################################################

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

			#############################
			# 11  Seq_Exp_ID ### SKIP ###
			#############################

			# 12  Genome_Ref=the reference genome used in the analysis pipeline

				REF_GENOME=${PLATFORM_UNIT_ARRAY[10]}

			###########################
			# 13  Operator ### SKIP ###
			##########################################
			# 14  Extra_VCF_Filter_Params ### SKIP ###
			##########################################

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

				FAMILY=${PLATFORM_UNIT_ARRAY[17]}

			# 22 DAD

				FAMILY=${PLATFORM_UNIT_ARRAY[17]}

			# 23 GENDER

				FAMILY=${PLATFORM_UNIT_ARRAY[17]}

			# 24 PHENOTYPE

				FAMILY=${PLATFORM_UNIT_ARRAY[17]}
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
			-S /usr/bin/env bash \
			-cwd \
			-V \
			-q $QUEUE_LIST \
			-p $PRIORITY \
		-N A.01-BWA"_"$SGE_SM_TAG"_"$FCID"_"$LANE"_"$INDEX \
			-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"_"$FCID"_"$LANE"_"$INDEX"-BWA.log" \
			-j y \
		$SCRIPT_DIR/A.01_BWA.sh \
			$BWA_DIR \
			$SAMBLASTER_DIR \
			$JAVA_1_8 \
			$PICARD_DIR \
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

###############################################################################
# create a hold job id qsub command line based on the number of ###############
# submit merging the bam files created by bwa mem above #######################
# only launch when every lane for a sample is done being processed by bwa mem #
# I want to clean this up eventually, but not in the mood for it right now. ###
###############################################################################
	#########################################################################################
	# I am setting the heap space and garbage collector threads now #########################
	# doing this does drastically decrease the load average ( the gc thread specification ) #
	#########################################################################################

		awk 1 $SAMPLE_SHEET \
			| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d; /^,/d' \
			| awk 'BEGIN {FS=","; OFS="\t"} NR>1 {print $1,$8,$2"_"$3"_"$4,$2"_"$3"_"$4".bam",$8,$10}' \
			| awk 'BEGIN {OFS="\t"} {sub(/@/,"_",$5)} {print $1,$2,$3,$4,$5,$6}' \
			| sort -k 1,1 -k 2,2 -k 3,3 -k 6,6 \
			| uniq \
			| $DATAMASH_DIR/datamash -s -g 1,2 collapse 3 collapse 4 unique 5 unique 6 \
			| awk 'BEGIN {FS="\t"} \
				gsub(/,/,",A.01-BWA_"$5"_",$3) \
				gsub(/,/,",INPUT=" "'$CORE_PATH'" "/" $1"/TEMP/",$4) \
				{print "qsub",\
				"-S /bin/bash",\
				"-cwd",\
				"-V",\
				"-q","'$QUEUE_LIST'",\
				"-p","'$PRIORITY'",\
				"-N","C.01-MARK_DUPLICATES_"$5"_"$1,\
				"-o","'$CORE_PATH'/"$1"/LOGS/"$2"/"$2"-MARK_DUPLICATES.log",\
				"-j y",\
				"-hold_jid","A.01-BWA_"$5"_"$3, \
				"'$SCRIPT_DIR'""/C.01_MARK_DUPLICATES.sh",\
				"'$JAVA_1_8'",\
				"'$PICARD_DIR'",\
				"'$SAMBAMBA_DIR'",\
				"'$CORE_PATH'",\
				$1,\
				$2,\
				"'$SAMPLE_SHEET'",\
				"'$SUBMIT_STAMP'",\
				$6,\
				"INPUT=" "'$CORE_PATH'" "/" $1"/TEMP/"$4"\n""sleep 0.1s"}'

	################################################################################
	# create an array at the SM tag level to populate aggregated sample variables. #
	################################################################################

		CREATE_SAMPLE_ARRAY ()
		{
			SAMPLE_ARRAY=(`awk 1 $SAMPLE_SHEET \
				| sed 's/\r//g; /^$/d; /^[[:space:]]*$/d' \
				| awk 'BEGIN {FS=","} $8=="'$SM_TAG'" {split($19,INDEL,";"); print $1,$5,$6,$7,$8,$9,$10,$12,$15,$16,$17,$18,INDEL[1],INDEL[2]}' \
				| sort \
				| uniq`)

			#  1  Project=the Seq Proj folder name
			PROJECT=${SAMPLE_ARRAY[0]}

			###################################################################
			#  2 SKIP : FCID=flowcell that sample read group was performed on #
			###################################################################

			############################################################################
			#  3 SKIP : Lane=lane of flowcell that sample read group was performed on] #
			############################################################################

			################################
			#  4 SKIP Index=sample barcode #
			################################

			#  5  Platform=type of sequencing chemistry matching SAM specification
			PLATFORM=${SAMPLE_ARRAY[1]}

			#  6  Library_Name=library group of the sample read group, Used during Marking Duplicates to determine if molecules are to be considered as part of the same library or not
			LIBRARY=${SAMPLE_ARRAY[2]}

			#  7  Date=should be the run set up date to match the seq run folder name, but it has been arbitrarily populated
			RUN_DATE=${SAMPLE_ARRAY[3]}

			#  8  SM_Tag=sample ID
			SM_TAG=${SAMPLE_ARRAY[4]}
			SGE_SM_TAG=$(echo $SM_TAG | sed 's/@/_/g') # If there is an @ in the qsub or holdId name it breaks

			#  9  Center=the center/funding mechanism
			CENTER=${SAMPLE_ARRAY[5]}

			# 10  Description=Generally we use to denote the sequencer setting (e.g. rapid run)
			# “HiSeq-X”, “HiSeq-4000”, “HiSeq-2500”, “HiSeq-2000”, “NextSeq-500”, or “MiSeq”.
			SEQUENCER_MODEL=${SAMPLE_ARRAY[6]}

			#############################
			# 11  Seq_Exp_ID ### SKIP ###
			#############################

			# 12  Genome_Ref=the reference genome used in the analysis pipeline
			REF_GENOME=${SAMPLE_ARRAY[7]}

			###########################
			# 13  Operator ### SKIP ###
			##########################################
			# 14  Extra_VCF_Filter_Params ### SKIP ###
			##########################################

			# 15  TS_TV_BED_File=where ucsc coding exons overlap with bait and target bed files
			TITV_BED=${SAMPLE_ARRAY[8]}

			# 16  Baits_BED_File=a super bed file incorporating bait, target, padding and overlap with ucsc coding exons.
			# Used for limited where to run base quality score recalibration on where to create gvcf files.
			BAIT_BED=${SAMPLE_ARRAY[9]}

				# since the mendel changes capture products need a way to define a 4th bed file which is the union of the different captures used.
					if [[ $PROJECT = "M_Valle"* ]];
						then
							HC_BAIT_BED=${MERGED_MENDEL_BED_FILE}
					elif [[ $PROJECT = "H_Cutting"* ]];
						then
							HC_BAIT_BED=${MERGED_CUTTING_BED_FILE}
					else
						HC_BAIT_BED=${BAIT_BED}
					fi

			# 17  Targets_BED_File=bed file acquired from manufacturer of their targets.
			TARGET_BED=${SAMPLE_ARRAY[10]}

			# 18  KNOWN_SITES_VCF=used to annotate ID field in VCF file. masking in base call quality score recalibration.
			DBSNP=${SAMPLE_ARRAY[11]}

			# 19  KNOWN_INDEL_FILES=used for BQSR masking, sensitivity in local realignment.
			KNOWN_INDEL_1=${SAMPLE_ARRAY[12]}
			KNOWN_INDEL_2=${SAMPLE_ARRAY[13]}
		}

	#############################################
	## using data only in the baited intervals ##
	#############################################
	## REMINDER TO HANDLE THE NEW JAR FILE NAME #
	#############################################

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

		FIX_BED_FILES ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
			-N A.00-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-FIX_BED_FILES.log" \
				-j y \
			-hold_jid C.01-MARK_DUPLICATES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/A.00_FIX_BED.sh \
				$SAMTOOLS_DIR \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$BAIT_BED \
				$TARGET_BED \
				$TITV_BED \
				$REF_GENOME
		}

	# run bqsr on the using bait bed file

		RUN_BQSR ()
		{
			echo \
			qsub \
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
			-N D.01-PERFORM_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-PERFORM_BQSR.log" \
				-j y \
			-hold_jid C.01-MARK_DUPLICATES"_"$SGE_SM_TAG"_"$PROJECT,A.00-FIX_BED_FILES"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/D.01_PERFORM_BQSR.sh \
				$JAVA_1_8 \
				$GATK_DIR_4011 \
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
				-S /bin/bash \
				-cwd \
				-V \
				-q $QUEUE_LIST \
				-p $PRIORITY \
			-N E.01-APPLY_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
				-o $CORE_PATH/$PROJECT/LOGS/$SM_TAG/$SM_TAG"-APPLY_BQSR.log" \
				-j y \
			-hold_jid D.01-PERFORM_BQSR"_"$SGE_SM_TAG"_"$PROJECT \
			$SCRIPT_DIR/E.01_APPLY_BQSR.sh \
				$JAVA_1_8 \
				$GATK_DIR_4011 \
				$CORE_PATH \
				$PROJECT \
				$SM_TAG \
				$REF_GENOME \
				$SAMPLE_SHEET \
				$SUBMIT_STAMP
		}

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
