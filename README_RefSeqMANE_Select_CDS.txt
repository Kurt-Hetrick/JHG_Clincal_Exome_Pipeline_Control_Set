16 June 2020: KNH

#######################################
########## REFSEQ SELECT CDS ##########
#######################################

	# 1. download file from (not that this is from the current directory, files will change and this will end up in a versioned directory)

			https://ftp.ncbi.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz

			Name                                         Last modified      Size  
			Parent Directory                                                  -   
			GRCh37_latest_assembly_report.txt            2019-10-24 18:09   28K  
			GRCh37_latest_clinvar.vcf.gz                 2020-06-16 12:06   28M  
			GRCh37_latest_dbSNP_all.vcf.gz               2020-05-14 13:10   15G  
			GRCh37_latest_genomic.fna.gz                 2019-10-24 18:09  900M  
			GRCh37_latest_genomic.gff.gz                 2019-10-24 18:10   25M  
			GRCh37_latest_genomic.gtf.gz                 2019-10-24 18:10   18M  
			GRCh37_latest_knownrefseq_alignments.bam     2019-10-24 18:11   49M  
			GRCh37_latest_knownrefseq_alignments.bam.bai 2019-10-24 18:11  2.0M  
			GRCh37_latest_protein.faa.gz                 2019-10-24 18:10   13M  
			GRCh37_latest_protein.gpff.gz                2019-10-24 18:10  109M  
			GRCh37_latest_rna.fna.gz                     2019-10-24 18:10   52M  
			GRCh37_latest_rna.gbff.gz                    2019-10-24 18:10  231M  
			README.txt                                   2019-11-05 15:49  2.2K  

	# 2. extract gene names and transcript ids for all CDS feature
	# there are 19212 refseq select transcripts that have cds.
	# 19197 are on the primary assembly
	# 9716 are forward strand
	# 9481 are reverse strand

		zegrep "RefSeq Select|MANE Select" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' >| refseq_select_cds_genes.txt

		zegrep "RefSeq Select|MANE Select" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /Parent=rna-(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' >| refseq_select_cds_transcripts.txt

	# 3. exon number forward strand transcripts. remove unplaced/unlocalized contigs.

		paste refseq_select_cds_genes.txt refseq_select_cds_transcripts.txt | awk 'BEGIN {OFS="\t"} {split($1,chrom,"."); print chrom[1],$2,$3,$5,$10,$4}' | sed 's/^NC_[0]*//g' | grep -v ^N | awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' | awk '$6=="+" {print $0 "\t" ++count[$5]}' | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' >| RefSeqSelect_annotated_transcripts_forward.txt

	# 4. exon number reverse strand transcripts. remove unplaced/unlocalized contigs. concatenate to forward strand transcripts.

		paste refseq_select_cds_genes.txt refseq_select_cds_transcripts.txt | awk 'BEGIN {OFS="\t"} {split($1,chrom,"."); print chrom[1],$2,$3,$5,$10,$4}' | sed 's/^NC_[0]*//g' | grep -v ^N | awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' | awk '$6=="-"' | sort -k4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr | awk '{print $0 "\t" ++count[$5]}' | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' | cat RefSeqSelect_annotated_transcripts_forward.txt - >| RefSeqSelect_annotated_transcripts_unsorted.txt

	# 5. reference genome sort annotated refseq select transcript file.
	# since some genes overlap, if you don't sort on gene name you will get a higher gene count.

		(awk '$1~/^[0-9]/' RefSeqSelect_annotated_transcripts_unsorted.txt | sort -k1,1n -k2,2n ; \
		awk '$1=="X"' RefSeqSelect_annotated_transcripts_unsorted.txt | sort -k 2,2n ; \
		awk '$1=="Y"' RefSeqSelect_annotated_transcripts_unsorted.txt | sort -k 2,2n) \
		>| GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed

##############################################################
########## WHAT IS IN OMIM NOT IN REFSEQ SELECT CDS ##########
##############################################################

	# Grab the uniq gene symbols from mim2gene

		grep -v ^# mim2gene_200611mbs.txt | awk 'BEGIN {FS="\t"} {print $4}' | sort | uniq | awk 'NR>1 {print  $1 "\t" "OMIM"}' >| mim2gene_200611mbs_uniq_gene_symbols.txt

	# Check to see if there are any genes in omim that are not in refseq select by gene symbol

		cut -f 4 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | sort | uniq | awk '{print $1 "\t" "REFSEQ"}' | cat mim2gene_200611mbs_uniq_gene_symbols.txt - | sort -k 1,1 | datamash -g 1 collapse 2 | awk '$2=="OMIM"' | wc -l
		928

		# there 930 rows in the final spreadsheet b/c IGH has 3 rows due to it having 3 OMIM ID #s
		# it's not kept

	## Create a text file for these.

		cut -f 4 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | sort | uniq | awk '{print $1 "\t" "REFSEQ"}' | cat mim2gene_200611mbs_uniq_gene_symbols.txt - | sort -k 1,1 | datamash -g 1 collapse 2 | awk '$2=="OMIM"' >| mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt

	# Sanity check the others (how many are just in refseq, how many are in omim and refseq)
		# 1. how many are only in refseq

			cut -f 4 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | sort | uniq | awk '{print $1 "\t" "REFSEQ"}' | cat mim2gene_200611mbs_uniq_gene_symbols.txt - | sort -k 1,1 | datamash -g 1 collapse 2 | awk '$2=="REFSEQ"' | wc -l
			3918

		# 2. how many are in are in both refseq and omim (check both ways)

			cut -f 4 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | sort | uniq | awk '{print $1 "\t" "REFSEQ"}' | cat mim2gene_200611mbs_uniq_gene_symbols.txt - | sort -k 1,1 | datamash -g 1 collapse 2 | awk '$2=="REFSEQ,OMIM"' | wc -l
			0

			cut -f 4 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | sort | uniq | awk '{print $1 "\t" "REFSEQ"}' | cat mim2gene_200611mbs_uniq_gene_symbols.txt - | sort -k 1,1 | datamash -g 1 collapse 2 | awk '$2=="OMIM,REFSEQ"' | wc -l
			15279

	# for the genes that are not in refseq select cds, grab them back out of the mim2gene file.
	# this is start of the excel spreadsheet for molly to review.
	# this has 930 rows in it b/c IGH is in there three times with different OMIM id #s.
	# it is not kept.

		for gene in $(cut -f 1 mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt) ; do grep -v ^# mim2gene_200611mbs.txt | awk 'BEGIN {FS="\t"} $4=="'$gene'" {print $0}' ; done >| OMIM_NOT_IN_REFSEQSELECT_CDS.txt

####################################################################################################
##### Annotate OMIM_NOT_IN_REFSEQSELECT_CDS.txt for molly to check what she wants ##################
####################################################################################################
##### add what gene symbols are present based on "Dbxref=GeneID:{#}," ##############################
##### what is on the cidr research twist bait and target bed files #################################
##### There are 63 genes that need to be found by Dbxref=GeneID b/c Entrez does not match HGNC #####
####################################################################################################

	# grab the Dbxref ID for omim genes that are not in refseq select cds.

		for gene in $(cut -f 1 mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt); do grep -w $gene mim2gene_200611mbs.txt | awk 'BEGIN {FS="\t";OFS="\t"} {print $4,$3}' ; done | sort | uniq >| mim2gene_200611mbs_uniq_gene_symbols_Entrez_IDs_notinRefSeqSelectCDS.txt

	# grab the entrez gene symbol for omim genes not in refseq select cds by their entrez Dbxref ID

		for gene in $(cut -f 2 mim2gene_200611mbs_uniq_gene_symbols_Entrez_IDs_notinRefSeqSelectCDS.txt) ; do zgrep "Dbxref=GeneID:$gene," /mnt/clinical/ddl/NGS/Exome_Resources/BED_FILES_TWIST/GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' | awk 'BEGIN {OFS="\t"} {split($1,chrom,"."); print chrom[1],$2,$3,$5,$10,$4}' | sed 's/^NC_[0]*//g' | grep -v ^N | awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}'; done >| mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS_ResSeqExon_byEntrezID.bed

	# How many omim genes not in RefSeq Select CDS are found by their Dbxref ID

		sort -k 1,1 -k 2,2n mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS_ResSeqExon_byEntrezID.bed | cut -f 4 | sort | uniq | wc -l
		760

	# How many omim genes not in RefSeq Select CDS are found by their Dbxref ID are in the cidr twist bait bed 

		sort -k 1,1 -k 2,2n mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS_ResSeqExon_byEntrezID_take2.bed | uniq | bedtools intersect -a - -b Baits_BED_File_TwistCUEXmito_20190415.bed | cut -f 4 | sort | uniq | wc -l
		371

	# How many omim genes not in RefSeq Select CDS are found by their Dbxref ID are in the cidr twist target bed 

		sort -k 1,1 -k 2,2n mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS_ResSeqExon_byEntrezID_take2.bed | uniq | bedtools intersect -a - -b Targets_BED_File_TwistCUEXmito_20190405.bed | cut -f 4 | sort | uniq | wc -l
		352

	# grab all the transcripts in refseq that are exon, cds, mrna or transcript features
	# that are omim genes that do not have a refseq select cds
	# only keep those that are on the primary assembly
	# didn't end up using this, but I think it was a good thing to do

		for gene in $(cut -f 2 mim2gene_200611mbs_uniq_gene_symbols_Entrez_IDs_notinRefSeqSelectCDS.txt) ; do zgrep "Dbxref=GeneID:$gene," /mnt/clinical/ddl/NGS/Exome_Resources/BED_FILES_TWIST/GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"||$3=="CDS"||$3=="mRNA"||$3=="transcript"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' | awk 'BEGIN {OFS="\t"} {split($1,chrom,"."); print chrom[1],$2,$3,$5,$10,$4}' | sed 's/^NC_[0]*//g' | grep -v ^N | awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}'; done >| all_the_features.txt

	# sanity check that i got them all

		sort -k 1,1 -k 2,2n all_the_features.txt | uniq | bedtools intersect -a - -b Baits_BED_File_TwistCUEXmito_20190415.bed | cut -f 4 | sort | uniq | wc -l
		372

		# the extra one is GGNBP1. truthfully, i don't know why I am missing this one above, but it's being dropped anyway b/c it is a potential pseudogene, so I'm not going to dig.

############################################################################################################
##### molly went through and highlighted those that she wanted to keep #####################################
##### I added a field of whether hgnc matched entrez #######################################################
##### made two text files. key files. ######################################################################
########## one for hgnc did not entrez and there was a refseq select transcript (34)########################
############### this will be used to swap entrez with hgnc symbols at the of this process ##################
########## another for hgnc did not match entrez and there was no ref seq select transcript (4) ############
############### this will be used to grab genes by exon feature and swap entrez with hgnc at the end #######
############################################################################################################
##### there are a total 335 genes that molly wants to keep #################################################
##### 3 are being dropped for various reasons (NOTCH2NLC,OFCC1,ATXN8). (332 left) ##########################
##### 34 hgnc do not match entrez, but have a refseq select transcript. changed at the end 298 left ########
##### 2 are not in the primary assembly for refseq grch37, but were in the old clinical bed file (ucsc) ####
############ those two (TEX28, KCNJ18) will be added at the end. 296 ################################
############ PSRR2 is a weird one, will be added. notes later. 295
##### that means that there are 295 genes to add ###########################################################
############ 11 have cds, hgnc matches entrez, but no refseq select transcript on the primary assembly #####
############ 276 do not have cds, hgnc matches entrez, but have exon #######################################
############ 4 are exon, no cds, no transcript
############ 4, hgnc does not have entrez and no cds #######################################################
############ there is one gene (ADSSL1) that does not transcript that ddl wants to keep (299) ##############
##### this makes a total of 19,496 genes (19197 refseq select cds transcripts) #############################
##### removing Y PAR leads to 19447 ########################################################################
############################################################################################################

##########################################################################
##### HGNC OMIM MATCHES ENTREZ. EXPLORE AND CHECK  #######################
##########################################################################

	# sanity check to make sure these don't have refseq select transcripts.
	# this is for genes, that molly wants to keep where the entrez ID matches the HGNC ID.
		# field 9 is HGNC_MATCHES_ENTREZ
		# field 11 is whether Molly wants to keep it or not.

			awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt | wc -l
			291

	# for omim genes where HGNC matches Entrez is there a CDS RefSeq select transcript.
	# answer: not on primary assembly

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="CDS"' ; done | egrep "RefSeq Select|MANE Select"

		# NT_167249.1     BestRefSeq      CDS     2357089 2357275 .       +       0       ID=cds-NP_001335178.1-2;Parent=rna-NM_001348249.1-2;Dbxref=GeneID:285834,Genbank:NP_001335178.1,HGNC:HGNC:27780,MIM:613918;Name=NP_001335178.1;Note=The RefSeq protein has 6 substitutions compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=CDS;gene=HCG22;inference=similar to AA sequence (same species):RefSeq:NP_001335178.1;product=protein PBMUCL2 precursor;protein_id=NP_001335178.1;tag=RefSeq Select
		# NT_167249.1     BestRefSeq      CDS     2358972 2359540 .       +       2       ID=cds-NP_001335178.1-2;Parent=rna-NM_001348249.1-2;Dbxref=GeneID:285834,Genbank:NP_001335178.1,HGNC:HGNC:27780,MIM:613918;Name=NP_001335178.1;Note=The RefSeq protein has 6 substitutions compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=CDS;gene=HCG22;inference=similar to AA sequence (same species):RefSeq:NP_001335178.1;product=protein PBMUCL2 precursor;protein_id=NP_001335178.1;tag=RefSeq Select

		# Molly Sheridan just FYI, doing some sanity checking, gene HCG22 does have a refseq select transcript. however the refseq select transcript is for an alternate haplotype of HLA. (ALT_REF_LOCI_7). For this, I have to pull out a non refseq select transcript on the primary assembly. - KNH

	# for when omim hgnc matches entrez gene are there any that cds.

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' | cut -f 5 ; done | uniq | wc -l
		12

	# for when omim hgnc matches entrez gene and have a cds, are they on the primary assembly.
	# HCG22 does not. it does have exon features on the primary assembly though.
	# so 11 genes where omim matches hgnc have cds.
	# so another fun fact is that out of these 11, 6 DO NOT have transcript IDs.
	# IGHD,IGHM,IGKC,TRAC,TRDC,TRGC1 DO NOT have transcript IDs
	# CCL3L1,GSTT1,KIR2DS4,KIR3DL1,SHANK3 do have transcript IDs

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' | cut -f 1,5 ; done | uniq
		NC_000017.10    CCL3L1
		NC_000022.10    GSTT1
		NT_167248.1     HCG22
		NT_167249.1     HCG22
		NC_000014.8     IGHD
		NW_004166863.1  IGHD
		NC_000014.8     IGHM
		NW_004166863.1  IGHM
		NC_000002.11    IGKC
		NC_000019.9     KIR2DS4
		NW_004166865.1  KIR2DS4
		NW_003571055.1  KIR2DS4
		NW_003571056.1  KIR2DS4
		NW_003571058.1  KIR2DS4
		NW_003571060.1  KIR2DS4
		NC_000019.9     KIR3DL1
		NW_004166865.1  KIR3DL1
		NW_003571060.1  KIR3DL1
		NC_000022.10    SHANK3
		NC_000014.8     TRAC
		NC_000014.8     TRDC
		NC_000007.13    TRGC1


#############################################################################################
########## ENTREZ MATCHES HGNC. NO REFSEQ SELECT. BUT HAVE CDS. HAVE TRANSCRIPT ID ##########
#############################################################################################

	# grab gene names of omim genes that match entrez gene names that have a cds but not a refseq select transcript. use transcripts on the primary assembly. hgnc matches entrez. select on CCL3L1,GSTT1,KIR2DS4,KIR3DL1,SHANK3 because these have transcript IDs

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | egrep "CCL3L1|GSTT1|KIR2DS4|KIR3DL1|SHANK3" | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,$4,foo[1]}' ; done >| OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_GENES.txt

	# grab transcripts name of omim genes that match entrez gene names that have a cds but not a refseq select transcript. use transcripts on the primary assembly. hgnc matches entrez.

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | egrep "CCL3L1|GSTT1|KIR2DS4|KIR3DL1|SHANK3" | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /Parent=rna-(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,$4,foo[1]}' ; done >| OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_TRANSCRIPTS.txt

	# when hgnc matches omim and there is a cds grab the forward strand genes.
	# filter to the transcript with the most defined bases
	# then with largest different in 3' and 5' ends
	# then to the one with "newest" accession id
	# number the exons for each transcript

		paste OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_GENES.txt OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_TRANSCRIPTS.txt \
		| sort -k 5,5 -k 10,10 -k 1,1 -k 2,2n -k 3,3n \
		| awk '{print $0 "\t" $3-$2}' \
		| datamash -g 1,5,10,4 first 2 last 3 sum 11 \
		| sort -k 2,2 -k 7,7nr -k 3,3r \
		| awk '{print $0 "\t" $6-$5}' \
		| sort -k 2,2 -k 7,7nr -k 8,8nr -k 3,3r \
		| datamash -g 1,2,4 first 3 \
		| for transcript in $(cut -f 4 -) ; \
		do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz \
		| awk 'BEGIN {FS="\t"} $3=="CDS"' \
		| awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' \
		| awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="+" {print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
		done \
		>| OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_FORWARD.txt

	# when hgnc matches omim and there is a cds grab the reverse strand genes.
	# filter to the transcript with the most defined bases
	# then with largest different in 3' and 5' ends
	# then to the one with "newest" accession id
	# number the exons for each transcript

		paste OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_GENES.txt OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_TRANSCRIPTS.txt \
		| sort -k 5,5 -k 10,10 -k 1,1 -k 2,2n -k 3,3n \
		| awk '{print $0 "\t" $3-$2}' \
		| datamash -g 1,5,10,4 first 2 last 3 sum 11 \
		| sort -k 2,2 -k 7,7nr -k 3,3r \
		| awk '{print $0 "\t" $6-$5}' \
		| sort -k 2,2 -k 7,7nr -k 8,8nr -k 3,3r \
		| datamash -g 1,2,4 first 3 \
		| for transcript in $(cut -f 4 -) ; \
		do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz \
		| awk 'BEGIN {FS="\t"} $3=="CDS"' \
		| awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' \
		| awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="-"' \
		| sort -k 4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr \
		| awk '{print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
		done \
		| cat OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_FORWARD.txt - \
		>| OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_COMPLETE_UNSORTED.txt

######################################################################################################
########## ENTREZ MATCHES HGNC. NO REFSEQ SELECT. BUT HAVE CDS. DOES NOT HAVE TRANSCRIPT ID ##########
######################################################################################################

	# grab gene names of omim genes that match entrez gene names that have a cds but not a refseq select transcript. use transcripts on the primary assembly. hgnc matches entrez. select on IGHD,IGHM,IGKC,TRAC,TRDC,TRGC1 because these DO NOT have transcript IDs
	# making the transcript name for these #N/A b/c that is what DDL did for ADSSL1
	# do the forward strand first

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; \
		do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz \
		| egrep "IGHD|IGHM|IGKC|TRAC|TRDC|TRGC1" \
		| awk 'BEGIN {FS="\t"} $3=="CDS"' \
		| awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' \
		| awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],"#N/A",$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="+" {print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
		done \
		>| OMIM_CDS_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_FORWARD.txt

	# now the reverse strand and combine with the forward strand

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; \
		do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz \
		| egrep "IGHD|IGHM|IGKC|TRAC|TRDC|TRGC1" \
		| awk 'BEGIN {FS="\t"} $3=="CDS"' \
		| awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' \
		| awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],"#N/A",$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="-"' \
		| sort -k 4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr \
		| awk '{print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
		done \
		| cat OMIM_CDS_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_FORWARD.txt - \
		>| OMIM_CDS_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_COMPLETE_UNSORTED.txt

################################################################################################
########## ENTREZ MATCHES HGNC. NO REFSEQ SELECT. NO CDS. USE EXON. HAS TRANSCRIPT ID ##########
################################################################################################

	# grab gene names of omim genes that match entrez gene names that DO NOT have a cds. use exon. use transcripts on the primary assembly. hgnc matches entrez. adding IGHA2|IGHE|IGLC1 to be excluded since these also do not have transcript IDs

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt | egrep -vw "KIR2DS4|KIR3DL1|SHANK3|CCL3L1|GSTT1|IGHD|IGHM|IGKC|TRAC|TRDC|TRGC1|IGHA2|IGHE|IGLC1") ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,$4,foo[1]}' ; done >| OMIM_EXON_NOT_REFSEQSELECT_GENES.txt

	# grab transcript names of omim genes that match entrez gene names that DO NOT have a cds. use exon. use transcripts on the primary assembly. hgnc matches entrez. adding IGHA2|IGHE|IGLC1 to be excluded since these also do not have transcript IDs

		for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt | egrep -vw "KIR2DS4|KIR3DL1|SHANK3|CCL3L1|GSTT1|IGHD|IGHM|IGKC|TRAC|TRDC|TRGC1|IGHA2|IGHE|IGLC1") ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"' | awk 'BEGIN {FS="\t"} match($9, /Parent=rna-(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,$4,foo[1]}' ; done > OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt

	# when hgnc matches omim and there is no cds, using exon. grab the forward strand genes. filter to the largest transcript and then the latest accession IDs. number the exons for each transcript.
	# double check below using transcript_id=
	# looks Parent=rna- is better because it looks like it would refer to the primary assembly. transcript_id can be the same b/w transcripts in alternate locations.

		# zgrep "NR_030321.1" GRCh37_latest_genomic.gff.gz | awk '$3=="exon"'

		# NC_000007.13    BestRefSeq      exon    73605528        73605624        .       +       .       ID=exon-NR_030321.1-1;Parent=rna-NR_030321.1;Dbxref=GeneID:693175,Genbank:NR_030321.1,HGNC:HGNC:32846,M
		# IM:615070,miRBase:MI0003602;gbkey=precursor_RNA;gene=MIR590;product=microRNA 590;transcript_id=NR_030321.1

		# NW_003871064.1  BestRefSeq      exon    1720434 1720530 .       +       .       ID=exon-NR_030321.1-2-1;Parent=rna-NR_030321.1-2;Dbxref=GeneID:693175,Genbank:NR_030321.1,HGNC:HGNC:32846,MIM:615070,mi
		# RBase:MI0003602;gbkey=precursor_RNA;gene=MIR590;product=microRNA 590;transcript_id=NR_030321.1

			paste OMIM_EXON_NOT_REFSEQSELECT_GENES.txt OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt \
			| sort -k 5,5 -k 10,10 -k 1,1 -k 2,2n -k 3,3n \
			| awk '{print $0 "\t" $3-$2}' \
			| datamash -g 1,5,10,4 first 2 last 3 sum 11 \
			| sort -k 2,2 -k 7,7nr -k 3,3r \
			| awk '{print $0 "\t" $6-$5}' \
			| sort -k 2,2 -k 7,7nr -k 8,8nr -k 3,3r \
			| datamash -g 1,2,4 first 3 \
			| for transcript in $(cut -f 4 -) ; \
			do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz \
			| awk 'BEGIN {FS="\t"} $3=="exon"' \
			| awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' \
			| awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' \
			| sed 's/^NC_[0]*//g' \
			| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
			| awk '$6=="+" {print $0 "\t" ++count[$5]}' \
			| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
			done \
			>| OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_FORWARD.txt

	# when hgnc matches omim and there is no cds, using exon. grab the reverse strand genes. filter to the largest transcript and then the latest accession IDs. number the exons for each transcript. combine with the forward strand genes formatted the same way.

		paste OMIM_EXON_NOT_REFSEQSELECT_GENES.txt OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt \
		| sort -k 5,5 -k 10,10 -k 1,1 -k 2,2n -k 3,3n \
		| awk '{print $0 "\t" $3-$2}' \
		| datamash -g 1,5,10,4 first 2 last 3 sum 11 \
		| sort -k 2,2 -k 7,7nr -k 3,3r \
		| awk '{print $0 "\t" $6-$5}' \
		| sort -k 2,2 -k 7,7nr -k 8,8nr -k 3,3r \
		| datamash -g 1,2,4 first 3 \
		| for transcript in $(cut -f 4 -) ; \
		do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz \
		| awk 'BEGIN {FS="\t"} $3=="exon"' \
		| awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' \
		| awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="-"' \
		| sort -k 4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr \
		| awk '{print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
		done \
		| cat OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_FORWARD.txt - \
		>| OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_COMPLETE_UNSORTED.txt

			# paste OMIM_EXON_NOT_REFSEQSELECT_GENES.txt OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt | datamash -g 1,5,10,4 first 2 last 3 | awk '{print $0 "\t" $6-$5}' | awk 'BEGIN {OFS="\t"} {if ($7<0) print $1,$2,$3,$4,$5,$6,$7*(-1); else print $0}' | sort -k 2,2 -k 7,7nr -k 3,3r | datamash -g 1,2,4 first 3 | for transcript in $(cut -f 4 -) ; do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz | awk '$3=="exon"' | awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' | awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' | sed 's/^NC_[0]*//g' | awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' | awk '$6=="-"' | sort -k4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr | awk '{print $0 "\t" ++count[$5]}' | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; done | cat OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_FORWARD.txt - > OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_COMPLETE_UNSORTED.txt

##########################################################################################################
########## ENTREZ MATCHES HGNC. NO REFSEQ SELECT. NO CDS. USE EXON. DOES NOT HAVE TRANSCRIPT ID ##########
##########################################################################################################

	for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; \
	do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz \
	| egrep "IGHA2|IGHE|IGLC1" \
	| awk 'BEGIN {FS="\t"} $3=="exon"' \
	| awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' \
	| awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],"#N/A",$4}' \
	| sed 's/^NC_[0]*//g' \
	| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
	| awk '$6=="+" {print $0 "\t" ++count[$5]}' \
	| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
	done \
	>| OMIM_EXON_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_FORWARD.txt

	for gene in $(awk 'BEGIN {FS="\t"} $9=="YES"&&$11~"Yes" {print $4}' OMIM_NOT_IN_REFSEQSELECT_CDS_ANNOTATED_MBS200720_KNH_FIX.txt) ; \
	do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz \
	| egrep "IGHA2|IGHE|IGLC1" \
	| awk 'BEGIN {FS="\t"} $3=="exon"' \
	| awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' \
	| awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],"#N/A",$4}' \
	| sed 's/^NC_[0]*//g' \
	| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
	| awk '$6=="-"' \
	| sort -k 4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr \
	| awk '{print $0 "\t" ++count[$5]}' \
	| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
	done \
	| cat OMIM_EXON_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_FORWARD.txt - \
	>| OMIM_EXON_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_COMPLETE_UNSORTED.txt

#################################################################
##### HGNC DOES NOT MATCH ENTREZ. THERE IS NO REFSEQ SELECT #####
#################################################################

	# FOR OMIM GENES THAT DO MATCH ENTREZ. do any of them have CDS. (answer is no)

		for gene in $(awk 'NR>1 {print $2}' HGNC_does_not_match_ENTREZ_NOT_RefSeqSelect.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="CDS"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' | cut -f 1,5 ; done | uniq

	# FOR OMIM GENES THAT DO MATCH ENTREZ. what do they look like for exon.
	# answer. and they are on the primary assembly.

		for gene in $(awk 'NR>1 {print $2}' HGNC_does_not_match_ENTREZ_NOT_RefSeqSelect.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} {split($5,foo,";"); print $1,$2,$3,$4,foo[1]}' | cut -f 1,5 ; done | uniq
		NC_000012.11    CLLU1OS
		NC_000009.11    DEC1
		NC_000007.13    SSPO
		NC_000006.11    TCP10

	# grab gene names of omim genes that DO NOT match entrez gene names AND DO NOT have a cds. use exon. use transcripts on the primary assembly.

		for gene in $(awk 'NR>1 {print $2}' HGNC_does_not_match_ENTREZ_NOT_RefSeqSelect.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"' | awk 'BEGIN {FS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,$4,foo[1]}' ; done >| HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_GENES.txt

	# grab transcript names of omim genes that DO NOT match entrez gene names AND DO NOT have a cds. use exon. use transcripts on the primary assembly.

		for gene in $(awk 'NR>1 {print $2}' HGNC_does_not_match_ENTREZ_NOT_RefSeqSelect.txt) ; do zgrep "gene=$gene;" GRCh37_latest_genomic.gff.gz | awk 'BEGIN {FS="\t"} $3=="exon"' | awk 'BEGIN {FS="\t"} match($9, /Parent=rna-(.*);/, a) {print $1,$4-1,$5,$7,a[1]}' | awk 'BEGIN {OFS="\t"} $1~/^NC/ {split($5,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,$4,foo[1]}' ; done >| HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt

	# when hgnc matches omim and there is no cds, using exon. grab the forward strand genes. filter to the largest transcript and then the latest accession IDs. number the exons for each transcript.

		paste HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_GENES.txt HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt \
		| sort -k 5,5 -k 10,10 -k 1,1 -k 2,2n -k 3,3n \
		| awk '{print $0 "\t" $3-$2}' \
		| datamash -g 1,5,10,4 first 2 last 3 sum 11 \
		| sort -k 2,2 -k 7,7nr -k 3,3r \
		| awk '{print $0 "\t" $6-$5}' \
		| sort -k 2,2 -k 7,7nr -k 8,8nr -k 3,3r \
		| datamash -g 1,2,4 first 3 \
		| for transcript in $(cut -f 4 -) ; \
		do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz \
		| awk 'BEGIN {FS="\t"} $3=="exon"' \
		| awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' \
		| awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print "Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="+" {print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; done \
		>| HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_FORWARD.txt

	# when hgnc matches omim and there is no cds, using exon. grab the reverse strand genes. filter to the largest transcript and then the latest accession IDs. number the exons for each transcript. combine with the forward strand genes formatted the same way.

		paste HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_GENES.txt HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS.txt \
		| sort -k 5,5 -k 10,10 -k 1,1 -k 2,2n -k 3,3n \
		| awk '{print $0 "\t" $3-$2}' \
		| datamash -g 1,5,10,4 first 2 last 3 sum 11 \
		| sort -k 2,2 -k 7,7nr -k 3,3r \
		| awk '{print $0 "\t" $6-$5}' \
		| sort -k 2,2 -k 7,7nr -k 8,8nr -k 3,3r \
		| datamash -g 1,2,4 first 3 \
		| for transcript in $(cut -f 4 -) ; \
		do zgrep "Parent=rna-$transcript;" GRCh37_latest_genomic.gff.gz \
		| awk 'BEGIN {FS="\t"}  $3=="exon"' \
		| awk 'BEGIN {FS="\t";OFS="\t"} match($9, /gene=(.*);/, a) {print $1,$4-1,$5,$7,"'$transcript'",a[1]}' \
		| awk 'BEGIN {OFS="\t"} {split($6,foo,";"); split($1,chrom,"."); print chrom[1],$2,$3,foo[1],$5,$4}' \
		| sed 's/^NC_[0]*//g' \
		| awk 'BEGIN {FS="\t";OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6; else if ($1=="24") print 		"Y",$2,$3,$4,$5,$6; else print $0}' \
		| awk '$6=="-"' \
		| sort -k 4,4 -k 5,5 -k 1,1 -k 2,2nr -k 3,3nr \
		| awk '{print $0 "\t" ++count[$5]}' \
		| awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$5,$7,$6}' ; \
		done \
		| cat HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_FORWARD.txt - \
		>| HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_COMPLETE_UNSORTED.txt

##################################################################

# combine all of the files

	cat GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed OMIM_CDS_NOT_REFSEQSELECT_W_TRANSCRIPT_IDS_COMPLETE_UNSORTED.txt OMIM_CDS_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_COMPLETE_UNSORTED.txt OMIM_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_COMPLETE_UNSORTED.txt OMIM_EXON_NOT_REFSEQSELECT_WO_TRANSCRIPT_IDS_COMPLETE_UNSORTED.txt HGNC_NOT_EQUAL_ENTREZ_EXON_NOT_REFSEQSELECT_TRANSCRIPTS_COMPLETE_UNSORTED.txt GENES_TO_ADD_MANUALLY.bed | awk 'BEGIN {OFS="\t"} {if ($1=="X") print "23",$2,$3,$4,$5,$6,$7; else if ($1=="Y") print "24",$2,$3,$4,$5,$6,$7; else print $0}' | sort -k 1,1n -k 2,2n -k 3,3n | awk 'BEGIN {OFS="\t"} {if ($1=="23") print "X",$2,$3,$4,$5,$6,$7; else if ($1=="24") print "Y",$2,$3,$4,$5,$6,$7; else print $0}' >| GRCh37_16June2020_RefSeqSelect_OMIM_CDS_exon_primary_assembly_annotated_transcripts.bed

# check against DDL.
# note that ADSSL1 doesn't have a transcript name. so it actually isn't caught by this query.

	awk 'BEGIN {OFS="\t"} {split($5,foo,"."); print foo[1],"REFSEQ"}' GRCh37_16June2020_RefSeqSelect_OMIM_CDS_exon_primary_assembly_annotated_transcripts.bed | sort -k 1,1 | uniq | cat RefSeq.Unique.GRCh37.FINAL.DDL_ADDED.Transcript_Only.bed - | sort -k 1,1 | datamash -g 1 collapse 2 | awk '$2=="DDL"' | cut -f 1 | awk '{print "grep -w","\x22"$1"\x22","RefSeq.Unique.GRCh37.FINAL.DDL_ADDED.bed"}' | bash
	14      105190607       105190799       ADSSL1  #N/A    1       +
	10      123276832       123276977       FGFR2   NM_000141       8       -
	2       110922096       110922264       NPHP1   NM_000272       8       -
	8       145049344       145049537       PLEC    NM_000445       2       -
	19      10906736        10906875        DNM2    NM_001005360    10      +
	6       43738443        43738983        VEGFA   NM_001025366    1       +
	2       220150708       220150719       DNAJB2  NM_001039550    10      +
	14      53569748        53569769        DDHD1   NM_001160147    3       -
	2       86444161        86444222        REEP1   NM_001164732    5       -
	10      112667437       112667594       BBIP1   NM_001195304    3       -
	2       202134233       202134328       CASP8   NM_001228       4       +
	X       107682587       107682601       COL4A6  NM_001847       1       -
	12      56093653        56093773        ITGA7   NM_002206       5       -
	9       35683155        35683238        TPM2    NM_003289       9       -
	9       35684728        35684804        TPM2    NM_003289       6       -
	11      615259  615318  IRF7    NM_004031       1       -
	13      28197011        28197387        POLR1D  NM_015972       3       +
	14      35182879        35182882        CFL2    NM_021914       1       -
	11      72003469        72004410        CLPB    NM_030813       17      -
	11      72083969        72084058        CLPB    NM_030813       5       -
	11      72114010        72114096        CLPB    NM_030813       3       -
	11      72145518        72145692        CLPB    NM_030813       1       -
	6       152466621       152466690       SYNE1   NM_033071       137     -
	6       152630965       152630969       SYNE1   NM_033071       89      -
	6       152832705       152832726       SYNE1   NM_033071       6       -
	14      76542939        76543034        IFT43   NM_052873       4       +
	3       129179696       129179849       IFT122  NM_052985       5       +
	22      20780030        20780031        SCARF2  NM_153334       11      -
	22      20783627        20783642        SCARF2  NM_153334       9       -

###############################################3

# grab genes in old ddl file that I want to keep.
# manually add version numbers to RefSeq.Unique.GRCh37.FINAL.30JANUARY2020_DDL_TRANSCRIPTS_TO_KEEP.bed

	egrep -w "FGFR2|ITGA7|NPHP1|PLEC|POLR1D|VEGFA|ADSSL1" ../BED_FILES/RefSeq.Unique.GRCh37.FINAL.30JANUARY2020.bed > RefSeq.Unique.GRCh37.FINAL.30JANUARY2020_DDL_TRANSCRIPTS_TO_KEEP.bed

# replace genes from old ddl bed file that I wan't to replace the entrez transcript with the ddl transcript

	egrep -w -v "FGFR2|ITGA7|NPHP1|PLEC|POLR1D|VEGFA" GRCh37_16June2020_RefSeqSelect_OMIM_CDS_exon_primary_assembly_annotated_transcripts.bed | cat - RefSeq.Unique.GRCh37.FINAL.30JANUARY2020_DDL_TRANSCRIPTS_TO_KEEP.bed >| GRCh37_16June2020_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_annotated_transcripts_raw.bed

# remove genes that are in the Y PAR.

	bedtools subtract -a GRCh37_16June2020_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_annotated_transcripts_raw.bed -b Y_PAR.bed >| GRCh37_16June2020_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_annotated_transcripts_NoYpar_unsorted.bed

# sort this bed file by reference genome order

	(awk '$1~/^[0-9]/' GRCh37_16June2020_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_annotated_transcripts_NoYpar_unsorted.bed | sort -k1,1n -k2,2n ; \
	awk '$1=="X"' GRCh37_16June2020_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_annotated_transcripts_NoYpar_unsorted.bed | sort -k 2,2n ; \
	awk '$1=="Y"' GRCh37_16June2020_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_annotated_transcripts_NoYpar_unsorted.bed | sort -k 2,2n) \
	>| GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated.bed

###########################################
#####  REPLACE ENTREZ NAMES WITH HGNC #####
###########################################

# MAKE A COPY OF THE COMPLETED BED FILE

	cp -vf GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated.bed GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed

# making sure that the refseq select transcripts where Entrez does not match HGNC are present in the bed file.
	wc -l HGNC_does_not_match_ENTREZ_RefSeqSelect.txt
	35 HGNC_does_not_match_ENTREZ_RefSeqSelect.txt

	for gene in $(awk 'NR>1 {print $2}' HGNC_does_not_match_ENTREZ_RefSeqSelect.txt) ; do awk '$4=="'$gene'"' GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed ; done | cut -f 4 | uniq | wc -l
	34

# replace entrez with hgnc for refseq select transcripts.

	awk 'NR>1 {print "sed -i" , "\x27" "s/" "\x5C" "<" $2 "\x5C" ">/" $1 "/g" "\x27" , "GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed"}' HGNC_does_not_match_ENTREZ_RefSeqSelect.txt | bash

# do again with non refseq select transcripts

	awk 'NR>1 {print "sed -i" , "\x27" "s/" "\x5C" "<" $2 "\x5C" ">/" $1 "/g" "\x27" , "GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed"}' HGNC_does_not_match_ENTREZ_NOT_RefSeqSelect.txt | bash

# sanity checks

	egrep -w "ADPRS|CEP20|CFAP91|CFAP94|CIBAR2|CYRIB|DNAI3|DUSP29|ELAPOR1|ELAPOR2|FAM86C1P|GFUS|IRAG2|KASH5|MACIR|MIR9-1HG|MYG1|PEDS1|POLR1F|POLR1G|POLR1H|TAMALIN|BPNT2|CDIN1|CEP43|CFAP251|CIBAR1|DNAAF6|DYNC2I1|DYNC2I2|KATNIP|LORICRIN|CT45A2|MORF4|CLLU1-AS1|DELEC1|SSPOP|TCP10L3" GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed | wc -l
	493

	egrep -w "ADPRS|CEP20|CFAP91|CFAP94|CIBAR2|CYRIB|DNAI3|DUSP29|ELAPOR1|ELAPOR2|FAM86C1P|GFUS|IRAG2|KASH5|MACIR|MIR9-1HG|MYG1|PEDS1|POLR1F|POLR1G|POLR1H|TAMALIN|BPNT2|CDIN1|CEP43|CFAP251|CIBAR1|DNAAF6|DYNC2I1|DYNC2I2|KATNIP|LORICRIN|CT45A2|MORF4|CLLU1-AS1|DELEC1|SSPOP|TCP10L3" GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed | cut -f 4 | uniq | wc -l
	38

egrep -w "ADPRHL2|FOPNL|MAATS1|CASC1|FAM92B|FAM49B|WDR63|DUPD1|KIAA1324|KIAA1324L|FAM86C1|TSTA3|LRMP|CCDC155|C5orf30|C1orf61|C12orf10|TMEM189|TWISTNB|CD3EAP|ZNRD1|GRASP|IMPAD1|C15orf41|FGFR1OP|WDR66|FAM92A|PIH1D3|WDR60|WDR34|KIAA0556|LOR|CXorf40A|MRVI1|CLLU1OS|DEC1|SSPO|TCP10" GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | wc -l
370

egrep -w "ADPRS|CEP20|CFAP91|CFAP94|CIBAR2|CYRIB|DNAI3|DUSP29|ELAPOR1|ELAPOR2|FAM86C1P|GFUS|IRAG2|KASH5|MACIR|MIR9-1HG|MYG1|PEDS1|POLR1F|POLR1G|POLR1H|TAMALIN|BPNT2|CDIN1|CEP43|CFAP251|CIBAR1|DNAAF6|DYNC2I1|DYNC2I2|KATNIP|LORICRIN|CT45A2|MORF4" GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed | wc -l
370

egrep -w "ADPRHL2|FOPNL|MAATS1|CASC1|FAM92B|FAM49B|WDR63|DUPD1|KIAA1324|KIAA1324L|FAM86C1|TSTA3|LRMP|CCDC155|C5orf30|C1orf61|C12orf10|TMEM189|TWISTNB|CD3EAP|ZNRD1|GRASP|IMPAD1|C15orf41|FGFR1OP|WDR66|FAM92A|PIH1D3|WDR60|WDR34|KIAA0556|LOR|CXorf40A|MRVI1|CLLU1OS|DEC1|SSPO|TCP10" GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | cut -f 4 | uniq | wc -l
34

# target territory

	awk 's+=($3-$2) {print s}' GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed | tail -n 1
	33755068


################## TO MAKE BAIT BED FILE #########################

cp GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_annotated_copy.bed \
GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_HGNC_annotated.bed

cut -f 1-3 GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_HGNC_annotated.bed \
>| GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_HGNC_short.bed

cat Baits_BED_File_TwistCUEXmito_20190415.bed GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_HGNC_short.bed | sort  -k 1,1 -k 2,2n -k 3,3n | bedtools merge -i - >| Baits_BED_File_TwistCUEXmito_20190415_GRCh37_RefSeqSelect_OMIM_DDL_CDS_exon_primary_assembly_NoYpar_Short_202006730.bed

#############################################################################################
########################################## FINISHED #########################################
#############################################################################################

	##########################################################################################
	##### CHECK WHAT WAS IN OMIM BUT NOT IN REFSEQ SELECT CDS TO OLD CLINICAL BED FILES ######
	##########################################################################################

		# for the genes that are not in refseq select cds, grab them out of the old coding bed file.

			for gene in $(cut -f 1 mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt) ; do grep -w $gene /mnt/clinical/ddl/NGS/Exome_Resources/BED_FILES/RefSeq.Unique.GRCh37.FINAL.30JANUARY2020.bed; done | cut -f 4 | uniq | wc -l
			40

			for gene in $(cut -f 1 mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt) ; do grep -w $gene /mnt/clinical/ddl/NGS/Exome_Resources/BED_FILES/RefSeq.Unique.GRCh37.FINAL.30JANUARY2020.bed; done | cut -f 4 | uniq > mim2gene_NotInRefSeqSelectCDS_butInOldCodingFile.txt

		# how many of them are NOT in the twist targeted capture bed file.

			for gene in $(cut -f 1 mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt) ; do grep -w $gene /mnt/clinical/ddl/NGS/Exome_Resources/BED_FILES/RefSeq.Unique.GRCh37.FINAL.30JANUARY2020.bed; done | bedtools subtract -a - -b Targets_BED_File_TwistCUEXmito_20190405.bed | wc -l
			72

		# how many of them are NOT in the twist bait capture bed file.

			for gene in $(cut -f 1 mim2gene_200611mbs_uniq_gene_symbols_notinRefSeqSelectCDS.txt) ; do grep -w $gene /mnt/clinical/ddl/NGS/Exome_Resources/BED_FILES/RefSeq.Unique.GRCh37.FINAL.30JANUARY2020.bed; done | bedtools subtract -a - -b Baits_BED_File_TwistCUEXmito_20190415.bed

			21      11097542        11097647        BAGE    NM_001187       2       -
			21      11098723        11098737        BAGE    NM_001187       1       -
			2       73927748        73928432        NAT8B   NM_016347       1       -
			1       148277378       148277551       NBPF14  NM_015383       42      -
			1       148320333       148320506       NBPF14  NM_015383       24      -
			10      51361722        51361910        PARG    NM_003631       4       -
			10      51362800        51363787        PARG    NM_003631       3       -
			10      118401625       118401793       PNLIPRP2        NM_005396       13      +
			9       35657748        35658015        RMRP    NR_003051       1       -
			2       122288457       122288583       RNU4ATAC        NR_023343       1       +
			X       153516034       153516762       TEX28   NM_001205201    2       -

##########################
###### MISCELLANEOUS #####
##########################

#################### DDL #############################

bedtools subtract -a RefSeq.Unique.GRCh37.FINAL.30JANUARY2020.bed -b RefSeq.Unique.GRCh37.FINAL.19Feb2018.bed > ../BED_FILES_TWIST/RefSeq.Unique.GRCh37.FINAL.DDL_ADDED.bed

awk 'BEGIN {OFS="\t"} {print $5,"DDL"}' RefSeq.Unique.GRCh37.FINAL.DDL_ADDED.bed | sort -k 1,1 | uniq >| RefSeq.Unique.GRCh37.FINAL.DDL_ADDED.Transcript_Only.bed

awk 'BEGIN {OFS="\t"} {split($5,foo,"."); print foo[1],"REFSEQ"}' GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | sort -k 1,1 | uniq >| GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.transcript_only.bed

######################################################
#################### HGNC ############################
######################################################

# making sure that the refseq select transcripts where Entrez does not match HGNC are present in the bed file.

wc -l HGNC_does_not_match_ENTREZ_RefSeqSelect.txt
33 HGNC_does_not_match_ENTREZ_RefSeqSelect.txt

khetric1@c6100-8> for gene in $(awk 'NR>1 {print $2}' HGNC_does_not_match_ENTREZ_RefSeqSelect.txt) ; do awk '$4=="'$gene'"' GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed ; done | cut -f 4 | uniq | wc -l
32

# replace entrez with hgnc

khetric1@c6100-8> cp GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts_copy.bed

khetric1@c6100-8> awk 'NR>1 {print "sed -i" , "\x27" "s/" "\x5C" "<" $2 "\x5C" ">/" $1 "/g" "\x27" , "GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts_copy.bed"}' HGNC_does_not_match_ENTREZ_RefSeqSelect.txt | bash

# sanity checks

khetric1@c6100-8> egrep -w "ADPRS|CEP20|CFAP91|CFAP94|CIBAR2|CYRIB|DNAI3|DUSP29|ELAPOR1|ELAPOR2|FAM86C1P|GFUS|IRAG2|KASH5|MACIR|MIR9-1HG|MYG1|PEDS1|POLR1F|POLR1G|POLR1H|TAMALIN|BPNT2|CDIN1|CEP43|CFAP251|CIBAR1|DNAAF6|DYNC2I1|DYNC2I2|KATNIP|LORICRIN" GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts_copy.bed | wc -l
349

khetric1@c6100-8> egrep -w "ADPRHL2|FOPNL|MAATS1|CASC1|FAM92B|FAM49B|WDR63|DUPD1|KIAA1324|KIAA1324L|FAM86C1|TSTA3|LRMP|CCDC155|C5orf30|C1orf61|C12orf10|TMEM189|TWISTNB|CD3EAP|ZNRD1|GRASP|IMPAD1|C15orf41|FGFR1OP|WDR66|FAM92A|PIH1D3|WDR60|WDR34|KIAA0556|LOR" GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed | wc -l
349

khetric1@c6100-8> wc -l GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts*
  191563 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.bed
  191563 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts_copy.bed
   19197 GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts.transcript_only.bed
  402323 total

mv GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts_copy.bed \
GRCh37_16June2020_RefSeqSelect_sorted_CDS_primary_assembly_annotated_transcripts_HGNC.bed
