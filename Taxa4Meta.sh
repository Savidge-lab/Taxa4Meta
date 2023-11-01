#!/bin/bash

usage(){
echo "
Taxa4Meta pipeline
Written by Qinglong Wu, Ph.D. from Baylor College of Medicine (Qinglong.Wu@bcm.edu)

Description: Perform human gut microbiome meta-analysis, i.e. 16S rRNA gene amplicon data generated from different sequencing strategies.
	     The core concept bases on the accurate taxamonic binning instead of OTU interpretation during 16S amplicon meta-analysis.

Usage: Taxa4Meta.sh --path /path-to-project-dir/ --threads [integer]
		    --fastq_truncee [integer] --minReadLen [integer] --maxReadLen [integer]
                    --fastq_stripleft [integer] --fastq_stripright [integer]
		    --BLCA-GenusConf [integer] --BLCA-SpeciesConf [integer]
                    --ReverseComplement [string] --PairedEnd [string] --ReferenceMode [string]

       (Note: all positional parameters must be provided in the same order as indicated above)

Parameters:
--path			Absolute path to the project directory. Arrange data: /ProjectDir/RawReadProcess/SampleDir/RawSequenceFiles
--threads		Number of threads to be specified. (must be integer)
--fastq_truncee		Truncate sequences from left until the total expected error is not higher than the specified value;
				usually 1 for Illumina short reads; could be 2 for 454 longer reads. (must be integer)
--minReadLen            During quality filtering (after paired-end merging), discard sequenced reads if they are shorter than
				the specified length. (must be integer)
--maxReadLen            During quality filtering (after paired-end merging), truncate sequence reads to the specified length,
				and keep sequence reads shorter than the length specified. (must be integer)
--fastq_stripleft	The length of any other technical sequence plus primer sequence to re removed from left.  (must be integer)
--fastq_stripright	The length of any other technical sequence plus primer sequence to be removed from right. (must be integer)
--BLCA-GenusConf	During BLCA taxonomic assignment, confidence cut-off for sequence annotation at genus level. (must be integer)
--BLCA-SpeciesConf	During BLCA taxonomic assignment, confidence cut-off for sequence annotation at species level. (must be integer)
--ReverseComplement	Perform reverse complement of filtered reads before clustering (reduce run time); must be TRUE or FALSE.
--PairedEnd		Switch to the process mode for paired-end reads, otherwise remain the process
				for single-end reads; must be TRUE or FALSE.
--ReferenceMode		Perform chimeras removal using reference-based mode, otherwise perform chimeras removal
				using both de novo mode and reference-based mode;  must be TRUE or FALSE.

Suggested key parameters:
V1V3:	forward orientation (minReadLen: 200, maxReadLen: 450; GenusConf: 90, SpeciesConf: 60)
	reverse orientation (minReadLen: 300, maxReadLen: 450; GenusConf: 90, SpeciesConf: 60)
V3V5:	forward orientation (minReadLen: 250, maxReadLen: 450; GenusConf: 85, SpeciesConf: 50)
	reverse orientation (minReadLen: 300, maxReadLen: 450; GenusConf: 85, SpeciesConf: 50)
V4:	forward orientation (minReadLen: 200, maxReadLen: 250; GenusConf: 70, SpeciesConf: 50)
	reverse orientation (minReadLen: 200, maxReadLen: 250; GenusConf: 70, SpeciesConf: 50)
V6V9:	forward orientation (minReadLen: 300, maxReadLen: 450; GenusConf: 90, SpeciesConf: 50)
	reverse orientation (minReadLen: 250, maxReadLen: 450; GenusConf: 90, SpeciesConf: 50)
"
}

if [ -z "$1" ] || [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--help" ]
then
	usage
	exit
fi


#assign variables with the input of positional arguments
ProgramDir=$(dirname $(readlink -f $0))	#get the path for the script and its dependencies
WorkDir=$2
THREADS=$4
maxEE=$6
MinLen=$8
MaxLen=${10}
TrimLeft=${12}
TrimRight=${14}
GenusConf=${16}
SpeciesConf=${18}
ReverseComplement=${20}
PairedEnd=${22}
ReferenceMode=${24}

cd $WorkDir


###################################################################################################################################################
#loop every sample directory

if [ $PairedEnd == "TRUE" ]
then	#process with paired-end reads from Illumina
	cd $WorkDir
	#loop every sample
	for SampleDir in $WorkDir/RawReadProcess/*
	do
		cd $SampleDir

		#for SRR files download from NCBI database: SRRxxxxx_1.fastq & SRRxxxxx_2.fastq in SRRxxxxx directory
		#"those variables might be changed if the input sequence is different"
	        ForwardSeq=`ls *_R1_*`
	        ReverseSeq=`ls *_R2_*`
		SampleName=`basename $SampleDir`

		#merge the paired-end reads
		$ProgramDir/scripts/vsearch --threads $THREADS \
				            --fastq_mergepairs $ForwardSeq \
					    --reverse $ReverseSeq \
					    --fastqout ''$SampleName'_merged.fq' \
					    --fastqout_notmerged_fwd ''$SampleName'_unmerged_forward.fq'

		#quality filtering and trimming, and strip forward and reverse primers by position
		$ProgramDir/scripts/vsearch --threads $THREADS \
				  	    --fastq_filter ''$SampleName'_merged.fq' \
					    --fastaout ''$SampleName'_merged_maxEE.fa' \
                        	            --fastq_stripleft $TrimLeft \
                                	    --fastq_stripright $TrimRight \
					    --fastq_truncee $maxEE \
					    --fastq_trunclen_keep $MaxLen \
					    --fastq_minlen $MinLen

		#apply for unmerged forward sequence, try to keep many reads
	        $ProgramDir/scripts/vsearch --threads $THREADS \
        	                            --fastq_filter ''$SampleName'_unmerged_forward.fq' \
                	                    --fastaout ''$SampleName'_unmerged_forward_maxEE.fa' \
					    --fastq_stripleft $TrimLeft \
					    --fastq_truncee $maxEE \
                                            --fastq_minlen $MinLen \
                                            --fastq_maxlen $MaxLen

		#concatenate clean reads of merged and forward unmerged reads
		cat ''$SampleName'_merged_maxEE.fa' ''$SampleName'_unmerged_forward_maxEE.fa' \
		  > ''$SampleName'_merged_filtered_noprimers.fa'

	        #perform reverse complement of the sequences if necessary
	        if [ $ReverseComplement == "TRUE" ]
	        then
        	        $ProgramDir/scripts/vsearch --threads $THREADS \
                	                            --fastx_revcomp ''$SampleName'_merged_filtered_noprimers.fa' \
                        	                    --fastaout ''$SampleName'_merged_clean.fa'
	        else
        	        #generate the same output for the next command
 		        cp ''$SampleName'_merged_filtered_noprimers.fa' ../''$SampleName'_merged_filtered_noprimers.fa'
			mv ../''$SampleName'_merged_filtered_noprimers.fa' ''$SampleName'_merged_clean.fa'
	        fi

		#QIIME format: SampleID_SeqID such as PatientID-Day-Collection_SeqID or PatientID.Day.Collection_SeqID
		#UPARSE/VSEARCH format: SampleID.SeqID such as PatientID_Day_Collection.SeqID
		prefix=`echo $SampleName | sed 's/-/_/g'`

		#Dereplicate sequences at sample level to reduce run time for OTU clustering
		$ProgramDir/scripts/vsearch --threads $THREADS \
					    --derep_fulllength ''$SampleName'_merged_clean.fa' \
					    --output ''$SampleName'_merged_clean_derep.fa' \
					    --relabel ''$prefix'.' \
					    --strand plus --fasta_width 0 --sizeout
	done

else	#process with single-end reads from 454 or Illumina (could be merged paired-end reads)
	cd $WorkDir
	for SampleDir in $WorkDir/RawReadProcess/*
	do
        	cd $SampleDir

        	#"those variables might be changed if the input sequence is different"
        	SingleSeq=`ls *`        #VSEARCH accepts .fq.gz or .fastq.gz files
        	SampleName=`basename $SampleDir`

	        #quality filtering and trimming, and strip primers positionaly (pay attention to the order of parameters)
	        $ProgramDir/scripts/vsearch --threads $THREADS \
        	                            --fastq_filter $SingleSeq \
                	                    --fastaout ''$SampleName'_maxEE.fa' \
                        	            --fastq_stripleft $TrimLeft \
					    --fastq_stripright $TrimRight \
 	                           	    --fastq_truncee $maxEE \
					    --fastq_trunclen_keep $MaxLen \
					    --fastq_minlen $MinLen

	        #perform reverse complement of the sequences: common in 454 sequencing if necessary
	        if [ $ReverseComplement == "TRUE" ]
	        then
	                $ProgramDir/scripts/vsearch --threads $THREADS \
        	                                    --fastx_revcomp ''$SampleName'_maxEE.fa' \
                	                            --fastaout ''$SampleName'_maxEE_clean.fa'
	        else
	                #generate the same output for the next command
	                cp ''$SampleName'_maxEE.fa' ../''$SampleName'_maxEE.fa'
			mv ../''$SampleName'_maxEE.fa' ''$SampleName'_maxEE_clean.fa'
	        fi

	        #QIIME format: SampleID_SeqID such as PatientID-Day-Collection_SeqID or PatientID.Day.Collection_SeqID
	        #UPARSE/VSEARCH format: SampleID.SeqID such as PatientID_Day_Collection.SeqID
	        prefix=`echo $SampleName | sed 's/\./_/g' | sed 's/-/_/g'`

	        #Dereplicate sequences at sample level to reduce run time for OTU clustering
	        $ProgramDir/scripts/vsearch --threads $THREADS \
        	                            --derep_fulllength ''$SampleName'_maxEE_clean.fa' \
                	                    --output ''$SampleName'_maxEE_clean_derep.fa' \
                        	            --relabel ''$prefix'.' \
                                	    --strand plus --fasta_width 0 --sizeout
	done

fi


###################################################################################################################################################
#prepare for OTU clustering

mkdir $WorkDir/Taxa4Meta_output
mkdir $WorkDir/Taxa4Meta_output/temp && cd $WorkDir/Taxa4Meta_output/temp

if [ $PairedEnd == "TRUE" ]
then	#prcess if paired-end reads
	cat $WorkDir/RawReadProcess/*/*_merged_clean_derep.fa > all_maxEE_clean_derep.fa
	rm $WorkDir/RawReadProcess/*/*merged*	#save storage
else	#prcess if single-end reads
	cat $WorkDir/RawReadProcess/*/*_maxEE_clean_derep.fa > all_maxEE_clean_derep.fa
	rm $WorkDir/RawReadProcess/*/*maxEE*	#save storage
fi

#First dereplicate across samples and remove singletons
$ProgramDir/scripts/vsearch --threads $THREADS \
		            --sizein --sizeout --fasta_width 0 \
                            --derep_fulllength all_maxEE_clean_derep.fa \
			    --uc all.derep.uc \
			    --output all.derep.fasta

#Remove chimeras
if [ $ReferenceMode == "TRUE" ]
then
	#only reference-based chimera detection (SILVA 132)
	$ProgramDir/scripts/vsearch --threads $THREADS \
	                            --sizein --sizeout --fasta_width 0 \
        	                    --uchime_ref all.derep.fasta \
                	            --db $ProgramDir/databases/silva_132_99_16S.fna \
	                       	    --nonchimeras all.derep.nonchimeras.fasta
else
	#de novo chimera detection, time-consuming for large dataset
	$ProgramDir/scripts/vsearch --threads $THREADS \
	                            --sizein --sizeout --fasta_width 0 \
				    --uchime_denovo all.derep.fasta \
				    --nonchimeras all.denovo.nonchimeras.fasta

	#continue with reference-based chimera detection (SILVA 132)
	$ProgramDir/scripts/vsearch --threads $THREADS \
                                    --sizein --sizeout --fasta_width 0 \
                                    --uchime_ref all.denovo.nonchimeras.fasta \
                                    --db $ProgramDir/databases/silva_132_99_16S.fna \
                                    --nonchimeras all.derep.nonchimeras.fasta
	rm all.denovo.nonchimeras.fasta
fi

#Extract all non-chimeric seqs
perl $ProgramDir/scripts/map.pl all_maxEE_clean_derep.fa all.derep.uc all.derep.nonchimeras.fasta > all_maxEE_clean_derep_nochimeras.fa

rm all.derep.uc all.derep.fasta all.denovo.nonchimeras.fasta all.derep.nonchimeras.fasta

#VSEARCH will internally sort the sequences by decreasing the length, cluster at 97% or 100% and relabel with OTU_n
#generate OTU table and OTU sequences (longest sequence, "no size (remove --sizeout)" in the identifiers - required by BLCA script)
#must use option "--cluster_fast" for input sequences with variable read length from the same region
#options "--centroids" + "--cluster_fast" will output the longest seed read as representative OTU sequence
$ProgramDir/scripts/vsearch --threads $THREADS --id 0.99 \
			    --strand plus --fasta_width 0 --sizein \
                            --cluster_fast all_maxEE_clean_derep_nochimeras.fa \
                            --relabel OTU_ \
                            --centroids all.otus.fasta \
                            --otutabout all.otutab.txt

###################################################################################################################################################
#assign taxonomy by BLCA with NCBI 16S Microbial DB (downloaded on July 25th, 2019)

cd $WorkDir/Taxa4Meta_output/temp
mkdir BLCA_taxonomy

#BLCA current version just use single thread for the script, split input into several new files of relatively even size
$ProgramDir/scripts/seqkit split --by-part 30 all.otus.fasta --out-dir Splits   #split sequences into 30 parts for downstream parallel processing

for SplitSeq in $PWD/Splits/*.fasta
do
	SplitSeqName=`basename $SplitSeq .fasta`
	#make sure you have dependencies installed (clustalo and blastn) and should be in your PATH
        python $ProgramDir/scripts/BLCA/2.blca_main.py -i $SplitSeq \
                                                       -r $ProgramDir/scripts/BLCA/db/16SMicrobial.ACC.taxonomy \
                                                       -q $ProgramDir/scripts/BLCA/db/16SMicrobial \
		 				       --cvrset 0.99 --iset 99 \
                                                       -o $PWD/BLCA_taxonomy/''$SplitSeqName'_BLCA_out.txt' &  #keep in background
done

wait  #wait the completion of all the background runs

cat ./BLCA_taxonomy/*_BLCA_out.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt
rm ./BLCA_taxonomy/*_BLCA_out.txt && rm -rf Splits  #remove the sequence splits

#parse the BLCA output file (classified OTU sequences)
sed -n '/Unclassified/!p' ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_classified.txt

sed 's/\t/;/g' ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_classified.txt > otu_taxonomy.txt      #fix for delimiter for parsing

#filter each taxa based on confidence score and keep all sorted gene IDs from LCA output for merging purpose
awk -v var=$GenusConf -v OFS="\t" -F ';' '{if ($3 >= var) print $1,$2; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_domain.txt
awk -v var=$GenusConf -v OFS="\t" -F ';' '{if ($5 >= var) print $1,$4; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_phylum.txt
awk -v var=$GenusConf -v OFS="\t" -F ';' '{if ($7 >= var) print $1,$6; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_class.txt
awk -v var=$GenusConf -v OFS="\t" -F ';' '{if ($9 >= var) print $1,$8; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_order.txt
awk -v var=$GenusConf -v OFS="\t" -F ';' '{if ($11 >= var) print $1,$10; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_family.txt
awk -v var=$GenusConf -v OFS="\t" -F ';' '{if ($13 >= var) print $1,$12; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_genus.txt
awk -v var=$SpeciesConf -v OFS="\t" -F ';' '{if ($15 >= var) print $1,$14; else print $1,"unclassified"}' otu_taxonomy.txt > taxa_species.txt

#join all ranks together
join -j 1 -a 1 -a 2 -t $'\t' -o auto taxa_domain.txt taxa_phylum.txt | \
               join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_class.txt | \
               join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_order.txt | \
              join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_family.txt | \
               join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_genus.txt | \
             join -j 1 -a 1 -a 2 -t $'\t' -o auto - taxa_species.txt > temp

cat temp | cut -d$'\t' -f 1 > temp1     #extract OTU ID column

cat temp | cut -d$'\t' -f 2-8 | sed 's/\t/;/g' > temp2  #extract taxonomy columns and convert to semi-colon delimiter

paste temp1 temp2 | sed 's/unclassified/Other/g' \
		  | sed 's/superkingdom://g' | sed 's/phylum://g' | sed 's/class://g' | sed 's/order://g' \
                  | sed 's/family://g' | sed 's/genus://g' | sed 's/species://g' \
		  | sed 's/Not Available/Other/g' \
		  > ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_classified_parsed.txt	#to be the same style with IDTAXA taxonomy file and taxonomic collapsing

rm taxa_*.txt temp* otu_taxonomy.txt

###################################################################################################################################################
#IDTAXA classifier for unclassified OTU sequences
#extract the unclassified OTU sequences from BLCA output file
sed -n '/Unclassified/p' ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment.txt | cut -d$'\t' -f 1 > ./BLCA_taxonomy/unclassified_OTU_IDs.txt

xargs samtools faidx all.otus.fasta < ./BLCA_taxonomy/unclassified_OTU_IDs.txt > all.otus_BLCA-unclassified.fasta

#run IDTAXA with the training set of RDP_v16-mod_March2018.RData (threshold: 70)
#input file - all.otus_BLCA-unclassified.fasta
#output file - all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_original.txt
Rscript --vanilla $ProgramDir/scripts/IDTAXA.R $WorkDir/Taxa4Meta_output/temp $ProgramDir/databases/RDP_v16-mod_March2018.RData

#parse the IDTAXA output file which does fit QIIME format of taxonomy file
cat all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_original.txt \
	| sed -n '/subclass/p' | sed -n '/suborder/p' | sed 's/"//g' | sed 's/\\//g' \
	| cut -d ';' -f 1,5,8,11,17,23,26 > all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_temp1.txt	#extract lines with subclass & suborder

cat all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_original.txt \
        | sed -n '/subclass/p' | sed -n '/suborder/!p' | sed 's/"//g' | sed 's/\\//g' \
        | cut -d ';' -f 1,5,8,11,17,20,23 > all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_temp2.txt	#extract lines with subclass only

cat all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_original.txt \
        | sed -n '/subclass/!p' | sed -n '/suborder/!p' \
	| sed -n '/unclassified_Root/!p' | sed -n '/Viruses/!p' | sed -n '/Eukaryota/!p' \
	| sed 1d | sed 's/"//g' | sed 's/\\//g' \
        | cut -d ';' -f 1,5,8,11,14,17,20 > all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_temp3.txt	#extract Bacteria and Archaea

cat all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_temp*.txt | sed 's/rootrank;//g' \
	> all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_parsed.txt	#clean taxonomy file a little bit
sed -i.backup 's/;unclassified_.*//g' all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_parsed.txt	#remove meaningless taxa

rm all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_temp*.txt all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_parsed.txt.backup	#save storage

###################################################################################################################################################
#make OTU table with taxonomy assignments from BLCA and IDTAXA

cd $WorkDir/Taxa4Meta_output/temp

#concatenate BLCA and IDTAXA taxonomy file
cat ./BLCA_taxonomy/otu_seqs_BLCA_tax_assignment_classified_parsed.txt all.otus_BLCA-unclassified_taxonomy_IDTAXA_0.7_RDP-V16-modified_parsed.txt \
	> all.otus_BLCA+IDTAXA_taxonomy.txt

cat all.otus_BLCA+IDTAXA_taxonomy.txt | cut -d$'\t' -f 1 > OTUs_list_bacteria_archaea.txt
sed -n '1p' all.otutab.txt > temp_1st-line.txt
grep -w -F -f OTUs_list_bacteria_archaea.txt all.otutab.txt > temp_matched.txt

cat temp_1st-line.txt temp_matched.txt > otu_table_bacteria_without_taxa.txt
rm temp* OTUs_list_bacteria_archaea.txt

rm -rf ./BLCA_taxonomy

biom convert -i otu_table_bacteria_without_taxa.txt -o otu_table_bacteria_without_taxa.biom --table-type="OTU table" --to-hdf5        #generate OTU table

biom add-metadata -i otu_table_bacteria_without_taxa.biom -o otu_table_bacteria_with_BLCA_IDTAXA.biom \
		  --observation-metadata-fp all.otus_BLCA+IDTAXA_taxonomy.txt \
		  --observation-header OTUID,taxonomy

biom summarize-table -i otu_table_bacteria_with_BLCA_IDTAXA.biom -o otu_table_bacteria_with_BLCA_IDTAXA_stats.txt       #get sequencing depth per sample

biom convert -i otu_table_bacteria_with_BLCA_IDTAXA.biom -o otu_table_bacteria_with_BLCA_IDTAXA.txt --to-tsv --header-key taxonomy

#extract OTU sequences based on the list in "otu_table_bacteria_with_BLCA_IDTAXA.txt"
cat otu_table_bacteria_with_BLCA_IDTAXA.txt | cut -d$'\t' -f 1 | sed 1,2d > temp
xargs samtools faidx all.otus.fasta < temp > final_otu_seqs.fasta
rm temp


#simplify the output directory
cd $WorkDir/Taxa4Meta_output/temp
cp otu_table_bacteria_with_BLCA_IDTAXA.biom ../
cp otu_table_bacteria_with_BLCA_IDTAXA.txt ../
cp final_otu_seqs.fasta ../

#rm -rf $WorkDir/Taxa4Meta_output/temp
