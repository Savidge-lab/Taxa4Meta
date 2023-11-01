## Taxa4Meta pipeline

## Description

Perform human gut microbiome meta-analysis, i.e. 16S rRNA gene amplicon data generated from different sequencing strategies. The core concept bases on the accurate taxamonic binning instead of OTU interpretation during gut meta-analysis.

## Dependencies for this pipeline
1) BLCA and NCBI 16S RefSeq database for species annotation of amplicons 
2) Samtools for fasta manipulation 
3) biom-format and h5py for OTU matrix conversion
4) DECIPHIER (R package): IdTaxa and its pre-built RDP databse (v16; curated by IdTaxa) for taxonomic annotation down to genus level
5) VSEARCH: sequence quality control and sequence clustering

## Usage

Do not run Taxa4Meta for combined datasets generated from different 16S variable regions, it has to be used for each dataset first and then merge all collpased species (L7) profiles for downstream meta-analysis.

Taxa4Meta.sh --path /path-to-project-dir/ --threads [integer] --fastq_truncee [integer] --minReadLen [integer] --maxReadLen [integer] --fastq_stripleft [integer] --fastq_stripright [integer] --BLCA-GenusConf [integer] --BLCA-SpeciesConf [integer] --ReverseComplement [string] --PairedEnd [string] --ReferenceMode [string]

(Note: all positional parameters must be provided in the same order as indicated above)

## Options

--path: Absolute path to the project directory. Arrange data: /ProjectDir/RawReadProcess/SampleDir/RawSequenceFiles

--threads: Number of threads to be specified. (must be integer)

--fastq_truncee: Truncate sequences from left until the total expected error is not higher than the specified value; usually 1 for Illumina short reads; could be 2 for 454 longer reads. (must be integer)

--minReadLen: During quality filtering (after paired-end merging), discard sequenced reads if they are shorter than the specified length. (must be integer)

--maxReadLen: During quality filtering (after paired-end merging), truncate sequence reads to the specified length, and keep sequence reads shorter than the length specified. (must be integer)

--fastq_stripleft: The length of any other technical sequence plus primer sequence to re removed from left.  (must be integer)

--fastq_stripright: The length of any other technical sequence plus primer sequence to be removed from right. (must be integer)

--BLCA-GenusConf: During BLCA taxonomic assignment, confidence cut-off for sequence annotation at genus level. (must be integer)

--BLCA-SpeciesConf: During BLCA taxonomic assignment, confidence cut-off for sequence annotation at species level. (must be integer)

--ReverseComplement: Perform reverse complement of filtered reads before clustering (reduce run time); must be TRUE or FALSE.

--PairedEnd: Switch to the process mode for paired-end reads, otherwise remain the process for single-end reads; must be TRUE or FALSE.

--ReferenceMode: Perform chimeras removal using reference-based mode, otherwise perform chimeras removal using both de novo mode and reference-based mode;  must be TRUE or FALSE.


## Suggested confidence and read length for key 16S variable regions

V1V3:   forward orientation (minReadLen: 200, maxReadLen: 450; GenusConf: 90, SpeciesConf: 60); reverse orientation (minReadLen: 300, maxReadLen: 450; GenusConf: 90, SpeciesConf: 60)
        
V3V5:   forward orientation (minReadLen: 250, maxReadLen: 450; GenusConf: 85, SpeciesConf: 50); reverse orientation (minReadLen: 300, maxReadLen: 450; GenusConf: 85, SpeciesConf: 50)
        
V4:     forward orientation (minReadLen: 200, maxReadLen: 250; GenusConf: 70, SpeciesConf: 50); reverse orientation (minReadLen: 200, maxReadLen: 250; GenusConf: 70, SpeciesConf: 50)

V6V9:   forward orientation (minReadLen: 300, maxReadLen: 450; GenusConf: 90, SpeciesConf: 50); reverse orientation (minReadLen: 250, maxReadLen: 450; GenusConf: 90, SpeciesConf: 50)



## Notes
1) to run your own analysis, you could manipulate lines 89-91 of Taxa4Meta.sh to fit the format of file names of fastq files and sample names 
2) For some cases, the lineage between genus and species could be different due to different taxonomic confidence set-up for these two ranks. To get consistent lineage for species/genus annotated by BLCA for downstream analysis, just replace the lineage for species/genus in the final taxonomy with NCBI taxonomic lineage (16SMicrobial.ACC.taxonomy_modification_records.xlsx under 'scripts' directory; text mapping based procedure) if you want to interpret results at collapsed species (L7) level
