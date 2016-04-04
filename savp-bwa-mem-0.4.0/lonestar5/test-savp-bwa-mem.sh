#!/bin/bash

#SBATCH -J test_savp-bwa-mem
#SBATCH -o test_savp-bwa-mem.o%j
#SBATCH -A iPlant-Collabs
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 24

#----
# Script for creating tar bundle of reference fasta and its index files for savp pipeline
# Creates index files
# Bundles everything up
# Note that format for samtools and picard has changed since last version
#----

module load bwa/0.7.12
module load picard-tools/1.141
module load samtools/1.3

## INPUTS
data_a="SRR2601691_1.fastq"
data_b="SRR2601691_2.fastq"
barcode="SRR2601691"   # make the user input this, required string
reference_bundle="e-coli-K-12.tar"  
cleanup=false

## PRE-PROCESSING INPUTS

# gunzip input files if needed
data_ext="${data_a##*.}"
if [ ${data_ext} = "gz" ]; then gunzip ${data_a}; data_a=${data_a%.*}; fi
data_ext="${data_b##*.}"
if [ ${data_ext} = "gz" ]; then gunzip ${data_b}; data_b=${data_b%.*}; fi

# untar reference bundle
tar -xvf ${reference_bundle}
ref=`echo "$reference_bundle" | sed 's/.tar//g'`
reference=${ref}.fa
dictionary=${ref}.dict

# number of cores to use per function  (40 tested as optimal on LS5.  Hyper-threading allows up to 48)
n=40

## MAIN FUNCTIONS


## If the reads are not paired, simply run samse, otherwise build the sam index for the mates and align them with sampe
if [ -z ${data_b} ]
    then bwa mem -t ${n} ${reference} ${data_a} | samtools view -bhS -@${n} - | samtools sort -@${n} - -o out_fasta_sorted.bam
else
    bwa mem -t ${n} ${reference} ${data_a} ${data_b} | samtools view -bhS -@${n} - | samtools sort -@${n} - -o out_fasta_sorted.bam
fi

java -Xmx32g -jar $TACC_PICARD_DIR/picard.jar AddOrReplaceReadGroups INPUT=out_fasta_sorted.bam OUTPUT=${barcode}.mem.bam RGLB=1 RGPL=illumina RGPU=all RGSM=${barcode} VALIDATION_STRINGENCY=SILENT

samtools index ${barcode}.mem.bam

# Clean up
if ${cleanup}; then echo "Cleaning up input and intermediate files"
	rm ${data_a}; rm ${data_b}										# remove inputs
	rm out_fasta_sorted.bam							# remove intermediate files
	rm ${reference_bundle}; rm ${reference}; rm ${dictionary} 		# remove reference files
	rm ${reference}.fai; rm ${reference}.amb; rm ${reference}.ann; 	# remove reference files
	rm ${reference}.bwt; rm ${reference}.pac; rm ${reference}.sa	# remove reference files
else
	echo "Skipping clean up"
fi
