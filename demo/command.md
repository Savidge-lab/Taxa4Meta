## analysis procedures

Step 1 - create directory "RawReadProcess" and split three V4 paired-end sequencing samples

Step 2 - Taxa4Meta command: /path-to/Taxa4Meta.sh --path /path-to/demo/RawReadProcess --threads 66 --fastq_truncee 1 --minReadLen 200 --maxReadLen 250 --fastq_stripleft 31 --fastq_stripright 32 --BLCA-GenusConf 70 --BLCA-SpeciesConf 50 --ReverseComplement FALSE --PairedEnd TRUE --ReferenceMode TRUE

## expected output: 
final_otu_seqs.fasta: OTU sequences
otu_table_bacteria_with_BLCA_IDTAXA.biom: biom format of OTU table
otu_table_bacteria_with_BLCA_IDTAXA.txt: tab-delimted format of OTU table
temp: a folder containing all intermediate files.
