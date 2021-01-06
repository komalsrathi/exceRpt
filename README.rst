.. |date| date::

*************************
exceRpt modified pipeline
*************************

:authors: Rob Kitchen
:editor: Komal S. Rathi
:contact: rathik@email.chop.edu
:organization: DBHi, CHOP
:status: This is "work in progress"
:date: |date|

Introduction
============

This is the modified exceRpt smallRNA pipeline.

exceRpt_smallRNA -- This choreographs the processing, filtering, and alignment of a single smallRNA-seq sample. 

mergePipelineRuns.R -- This script will take as input a directory containing 1 or more subdirectories or zipfiles containing output from the pipeline above. In this way, results from 1 or more smallRNA-seq samples can be combined, several QC plots are generated, and the read-counts are normalised ready for downstream analysis by clustering and/or differential expression.

Please see the exceRpt mainpage (https://rkitchen.github.io/exceRpt) for instructions as to how to use the software.


Installation on linux
=====================

Clone existing repo
-------------------

.. code-block:: bash

	# go to working directory
	cd /mnt/isilon/xing_lab/aspera/rathik

	# clone directory
	git clone https://github.com/rkitchen/exceRpt.git


Download databases
------------------

.. code-block:: bash

	# download databases
	mkdir exceRpt-db && cd exceRpt-db

	# get core database
	wget http://homes.gersteinlab.org/people/rrk24/exceRpt/exceRptDB_v4_CORE.tgz
	tar -xvf exceRptDB_v4_CORE.tgz
	cd DATABASE
	(folders: adapters, randomBits.dat, STAR_Parameters_Endogenous_smallRNA.in, STAR_Parameters_Exogenous.in, UniVec)

	# get pre-compiled genome and transcriptome indices
	wget http://org.gersteinlab.excerpt.s3-website-us-east-1.amazonaws.com/exceRptDB_v4_hg38_lowmem.tgz
	tar -xvf exceRptDB_v4_hg38_lowmem.tgz

	# get exogenous miRNAs and rRNAs
	wget http://org.gersteinlab.excerpt.s3-website-us-east-1.amazonaws.com/exceRptDB_v4_EXOmiRNArRNA.tgz
	tar -xvf exceRptDB_v4_EXOmiRNArRNA.tgz
	(folders: miRBase, NCBI_taxonomy_taxdump, ribosomeDatabase)

	# exogenous genomes
	wget http://org.gersteinlab.excerpt.s3-website-us-east-1.amazonaws.com/exceRptDB_v4_EXOGenomes.tgz
	tar -xvf exceRptDB_v4_EXOGenomes.tgz
	(folders: Genomes_BacteriaFungiMammalPlantProtistVirus)

	# set db path
	export DATABASE_PATH=/mnt/isilon/xing_lab/aspera/rathik/exceRpt-db/DATABASE


Download test dataset
---------------------

.. code-block:: bash

	# download test smRNA-seq data
	mkdir test && cd test
	wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX010/SRX010851/SRR026761/SRR026761.sra

	# create output directory
	mkdir output && cd output

Installation of dependencies using conda
----------------------------------------

.. code-block:: bash

	# install conda
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -p /mnt/isilon/cbmi/variome/rathik/tools/miniconda

	# create conda environment for excerpt dependencies
	conda create --prefix /mnt/isilon/xing_lab/aspera/rathik/excerpt_env

	# activate environment
	conda activate /mnt/isilon/xing_lab/aspera/rathik/excerpt_env

	# install individual dependencies 
	conda install -c bioconda fastx_toolkit
	conda install -c bioconda bowtie2
	conda install -c bioconda samtools
	conda install -c bioconda fastqc
	conda install -c bioconda sra-tools
	conda install -c bioconda star=2.4.2a
	conda install -c bioconda bbmap

Test smallRNA pipeline
----------------------

.. code-block:: bash

	# run makefile on test dataset
	nohup bash run-test.sh &

Output files
============

.. code-block:: bash

	# let's check test output files/directories
	$ ls output/ 
	SRR026761.qcResult	| Text file containing a variety of QC metrics for this sample
	SRR026761	| Directory containing the complete set of output files for this sample
	SRR026761.stats	| Text file containing a variety of alignment statistics for this sample
	SRR026761.log	| Text file containing normal logging information for this run
	SRR026761.err	| Text file containing error logging information for this run
	SRR026761_CORE_RESULTS_v5.0.0.tgz	| Archive containing only the most commonly used results for this sample

	# let's check the directory containing the complete set of output files
	$ ls SRR026761/

	# All compatible alignments against the transcriptome after invoking the library priority
	endogenousAlignments_Accepted.txt.gz

	# Contains the ID(s) of the RNA annotations indexed in the fifth column of the .txt.gz file above
	endogenousAlignments_Accepted.dict

	# Alignments (ungapped) to the endogenous genome
	endogenousAlignments_genome_Aligned.out.bam 
	endogenousAlignments_genome_Log.final.out 
	endogenousAlignments_genome_Log.out 
	endogenousAlignments_genome_Log.progress.out
	endogenousAlignments_genome__STARtmp

	# Summary of the alignment characteristics for genome-mapped reads
	endogenousAlignments_genome_Aligned.out.bam.CIGARstats.txt

	# Transcriptome alignments (ungapped) of reads mapped to the genome
	endogenousAlignments_genomeMapped_transcriptome_Aligned.out.bam 
	endogenousAlignments_genomeMapped_transcriptome_Log.final.out 
	endogenousAlignments_genomeMapped_transcriptome_Log.out.gz
	endogenousAlignments_genomeMapped_transcriptome_Log.progress.out
	endogenousAlignments_genomeMapped_transcriptome_SJ.out.tab
	endogenousAlignments_genomeMapped_transcriptome__STARtmp
	endogenousAlignments_genomeMapped_transcriptome_Unmapped.R1.fastq.gz

	# Contains read-depth across all gencode transcripts
	endogenousAlignments_genomeMapped_transcriptome_Aligned.out.sorted.bam.coverage.txt 

	# Transcriptome alignments (ungapped) of reads not mapped to the genome
	endogenousAlignments_genomeUnmapped_transcriptome_Aligned.out.bam 
	endogenousAlignments_genomeUnmapped_transcriptome_Log.final.out 
	endogenousAlignments_genomeUnmapped_transcriptome_Log.out.gz
	endogenousAlignments_genomeUnmapped_transcriptome_Log.progress.out
	endogenousAlignments_genomeUnmapped_transcriptome__STARtmp
	endogenousAlignments_genomeUnmapped_transcriptome_Unmapped.R1.fastq.gz

	# Alignments to the UniVec and rRNA sequences
	filteringAlignments_UniVec_and_rRNA_Aligned.out.bam 
	filteringAlignments_UniVec_and_rRNA_Log.final.out 
	filteringAlignments_UniVec_and_rRNA_Log.out 
	filteringAlignments_UniVec_and_rRNA_Log.progress.out
	filteringAlignments_UniVec_and_rRNA_SJ.out.tab
	filteringAlignments_UniVec_and_rRNA__STARtmp

	# Read counts of each annotated RNA using sense/antisense alignments
	readCounts_circRNA_antisense.txt
	readCounts_circRNA_sense.txt
	readCounts_gencode_antisense_geneLevel.txt
	readCounts_gencode_antisense.txt
	readCounts_gencode_sense_geneLevel.txt
	readCounts_gencode_sense.txt
	readCounts_miRNAmature_sense.txt
	readCounts_miRNAprecursor_antisense.txt
	readCounts_miRNAprecursor_sense.txt
	readCounts_piRNA_antisense.txt
	readCounts_piRNA_sense.txt
	readCounts_tRNA_antisense.txt
	readCounts_tRNA_sense.txt

	# Reads remaining after each QC / filtering / alignment step
	SRR026761.clipped.fastq.gz
	SRR026761.clipped.REMOVEDRepeatReads.fastq.gz
	SRR026761.clipped.trimmed.fastq.gz
	SRR026761.clipped.trimmed.filtered.fastq.gz
	SRR026761.clipped.trimmed.filtered.noUniVecOrRiboRNA.fastq.gz

	# FastQC output both before and after UniVec/rRNA contaminant removal
	SRR026761.clipped.trimmed.filtered_fastqc.html
	SRR026761.clipped.trimmed.filtered_fastqc.zip
	SRR026761.clipped.trimmed.filtered.noUniVecOrRiboRNA_fastqc.html
	SRR026761.clipped.trimmed.filtered.noUniVecOrRiboRNA_fastqc.zip

	# Counts of the number of reads of each length following adapter removal
	SRR026761.clipped.trimmed.filtered.readLengths.txt

	# unsure
	SRR026761.clipped.trimmed.filtered.tmp2
	SRR026761.clipped.trimmed.filtered.tmp.log
	SRR026761.clipped.trimmed.filtered.tmp.PhiX.log

	# Read counts mapped to UniVec & rRNA (and calibrator oligo, if used) sequences
	SRR026761.clipped.trimmed.filtered.UniVec_and_rRNA.counts
	SRR026761.clipped.trimmed.filtered.UniVec_and_rRNA.readCount 

	$ tar -tvf SRR026761_CORE_RESULTS_v5.0.0.tgz (same files as above just a collection of all important files)

	# Read counts of each annotated RNA using sense alignments: *sense.txt
	# Read counts of each annotated RNA using antisense alignments: *antisense.txt

	# Summary of the alignment characteristics for genome-mapped reads
	endogenousAlignments_genome_Aligned.out.bam.CIGARstats.txt
	# Contains read-depth across all gencode transcripts
	endogenousAlignments_genomeMapped_transcriptome_Aligned.out.sorted.bam.coverage.txt
	# circRNA
	readCounts_circRNA_antisense.txt
	readCounts_circRNA_sense.txt
	# gencode gene-level
	readCounts_gencode_antisense_geneLevel.txt
	readCounts_gencode_sense_geneLevel.txt
	# gencode
	readCounts_gencode_antisense.txt
	readCounts_gencode_sense.txt
	# miRNAmature
	readCounts_miRNAmature_sense.txt
	# miRNAprecursor
	readCounts_miRNAprecursor_antisense.txt
	readCounts_miRNAprecursor_sense.txt
	# piRNA
	readCounts_piRNA_antisense.txt
	readCounts_piRNA_sense.txt
	# tRNA
	readCounts_tRNA_antisense.txt
	readCounts_tRNA_sense.txt
	# FastQC output both before and after UniVec/rRNA contaminant removal
	SRR026761.clipped.trimmed.filtered_fastqc.zip
	SRR026761.clipped.trimmed.filtered.noUniVecOrRiboRNA_fastqc.zip
	# Counts of the number of reads of each length following adapter removal
	SRR026761.clipped.trimmed.filtered.readLengths.txt
	# Read counts mapped to UniVec & rRNA (and calibrator oligo, if used) sequences
	SRR026761.clipped.trimmed.filtered.UniVec_and_rRNA.counts
	# Text file containing normal logging information for this run
	SRR026761.log
	# Text file containing a variety of alignment statistics for this sample
	SRR026761.stats
	# Text file containing a variety of QC metrics for this sample
	SRR026761.qcResult

Note:
-----

Could not find the following files, may not be generated in the updated version:

[sampleID]/readCounts_miRNAmature_antisense.txt
[sampleID]/[sampleID].\*.knownAdapterSeq      | 3' adapter sequence guessed (from known adapters) in this sample
[sampleID]/[sampleID].\*.adapterSeq           | 3' adapter used to clip the reads in this run
[sampleID]/[sampleID].\*.qualityEncoding      | PHRED encoding guessed for the input sequence reads 


Script to merge pipeline runs
=============================

Install R packages
------------------

.. code-block:: bash

	conda install -c bioconda bioconductor-rgraphviz
	conda install -c r r-tidyverse
	conda install -c bioconda bioconductor-marray
	conda install -c conda-forge r-plyr
	conda install -c r r-scales
	conda install -c r r-reshape2
	conda install -c conda-forge r-gplots
	
Running the script
------------------

.. code-block:: bash

	Rscript mergePipelineRuns.R /mnt/isilon/xing_lab/aspera/rathik/urine_smallrna_output

Final output files:
-------------------

.. code-block:: bash

	exceRpt_adapterSequences.txt
	exceRpt_biotypeCounts.txt 
	exceRpt_circularRNA_ReadCounts.txt
	exceRpt_circularRNA_ReadsPerMillion.txt 
	exceRpt_DiagnosticPlots.pdf 
	exceRpt_exogenous_miRNA_ReadCounts.txt
	exceRpt_exogenous_miRNA_ReadsPerMillion.txt 
	exceRpt_exogenousRibosomal_taxonomyCumulative_ReadCounts.txt
	exceRpt_exogenousRibosomal_taxonomyCumulative_ReadsPerMillion.txt 
	exceRpt_exogenousRibosomal_taxonomySpecific_ReadCounts.txt
	exceRpt_exogenousRibosomal_taxonomySpecific_ReadsPerMillion.txt 
	exceRpt_exogenousRibosomal_TaxonomyTrees_aggregateSamples.pdf 
	exceRpt_exogenousRibosomal_TaxonomyTrees_perSample.pdf
	exceRpt_gencode_ReadCounts.txt
	exceRpt_gencode_ReadsPerMillion.txt 
	exceRpt_miRNA_ReadCounts.txt
	exceRpt_miRNA_ReadsPerMillion.txt 
	exceRpt_piRNA_ReadCounts.txt
	exceRpt_piRNA_ReadsPerMillion.txt 
	exceRpt_QCresults.txt 
	exceRpt_ReadLengths.txt 
	exceRpt_readMappingSummary.txt
	exceRpt_sampleGroupDefinitions.txt
	exceRpt_smallRNAQuants_ReadCounts.RData 
	exceRpt_smallRNAQuants_ReadsPerMillion.RData
	exceRpt_tRNA_ReadCounts.txt 
	exceRpt_tRNA_ReadsPerMillion.txt

