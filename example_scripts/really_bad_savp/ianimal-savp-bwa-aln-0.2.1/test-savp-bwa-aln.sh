#!/bin/bash

#SBATCH -J test_savp_bwa_aln
#SBATCH -o test_savp_bwa_aln.o%j
#SBATCH -A iPlant-Collabs
#SBATCH -p normal
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 16

module load picard
module load samtools
module load bwa

## INPUTS
data_a="HO.2081.AP.03.1_Unique.fastq"
data_b="HO.2081.AP.03.2_Unique.fastq"
barcode="HO.2081.AP.03"   # make the user input this, required string
reference_bundle="umd_3_1_Y_Mito.tar"  

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


## MAIN FUNCTIONS

## Generate index for forward reads
bwa aln -t 15 ${reference} ${data_a} > aln_a.sai

## If the reads are not paired, simpy run samse, otherwise build the sam index for the mates and align them with sampe
if [ -z ${data_b} ]
    then bwa samse ${reference} aln_a.sai ${data_a} | samtools view -bhS -@15 - | samtools sort -@15 - out_fasta_sorted
else
    bwa aln -t 15 ${reference} ${data_b} > aln_b.sai
    bwa sampe ${reference} aln_a.sai aln_b.sai ${data_a} ${data_b} | samtools view -bhS -@15 - | samtools sort -@15 - out_fasta_sorted
fi

java -Xmx12g -jar $TACC_PICARD_DIR/AddOrReplaceReadGroups.jar INPUT=out_fasta_sorted.bam OUTPUT=${barcode}.aln.bam RGLB=1 RGPL=illumina RGPU=all RGSM=${barcode} VALIDATION_STRINGENCY=SILENT

samtools index ${barcode}.aln.bam


