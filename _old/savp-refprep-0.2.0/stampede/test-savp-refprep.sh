#!/bin/bash

#SBATCH -J test_refprep
#SBATCH -o refprep_test.o%j
#SBATCH -A iPlant-Collabs
#SBATCH -p normal
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -n 16

#----
# Script for creating tar bundle of reference fasta and its index files for savp pipeline
# Creates index files
# Bundles everything up
#----

ref="foo.fasta.gz"

# is input compressed? if so, uncompress and update with new filename
ref_ext="${ref##*.}"
if [ ${ref_ext} = "gz" ]; then gunzip ${ref}; ref=${ref%.*}; fi

# standardize extension to .fa
rb="${ref%%.*}"   # reference basename
mv ${ref} ${rb}.fa
ref=${rb}.fa

gzip_output=1
if [ ${gzip_output} = 1 ]; then compress_flag="z"; else compress_flag=""; fi

module load bwa/0.7.7
module load picard
module load samtools

samtools faidx ${ref}
bwa index ${ref}
java -Xmx3g -jar $TACC_PICARD_DIR/CreateSequenceDictionary.jar R= ${ref} O= ${rb}.dict
tar -c${compress_flag}f ${rb}.tar ${ref} ${ref}.fai ${rb}.dict ${ref}.amb ${ref}.ann ${ref}.bwt ${ref}.pac ${ref}.sa

rm ${ref}; rm ${ref}.fai; rm ${rb}.dict; rm ${ref}.amb; rm ${ref}.ann; rm ${ref}.bwt; rm ${ref}.pac; rm ${ref}.sa
