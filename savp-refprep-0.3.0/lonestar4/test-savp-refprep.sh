#!/bin/bash
#$ -V
#$ -cwd # Start job in submission directory
#$ -N test-savp-refprep # Job Name
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12
#$ -q development       # Queue name normal
#$ -A iPlant-Collabs
#$ -l h_rt=01:00:00     # Run time (hh:mm:ss) 
#$ -M jcarson@tacc.utexas.edu   # Address for email notification
#$ -m be        # Email at Begin and End of job

#----
# Script for creating tar bundle of reference fasta and its index files for savp pipeline
# Creates index files
# Bundles everything up
#----

ref="foo.fasta.gz"
cleanup=true

# is input compressed? if so, uncompress and update with new filename
ref_ext="${ref##*.}"
if [ ${ref_ext} = "gz" ]; then gunzip ${ref}; ref=${ref%.*}; fi

# standardize extension to .fa
rb="${ref%%.*}"   # reference basename
mv ${ref} ${rb}.fa
ref=${rb}.fa

gzip_output=1
if [ ${gzip_output} = 1 ]; then compress_flag="z"; else compress_flag=""; fi

module load jdk64
module load bwa/0.7.7
module load picard
module load samtools

samtools faidx ${ref}
bwa index ${ref}
java -Xmx3g -jar $TACC_PICARD_DIR/CreateSequenceDictionary.jar R= ${ref} O= ${rb}.dict
tar -c${compress_flag}f ${rb}.tar ${ref} ${ref}.fai ${rb}.dict ${ref}.amb ${ref}.ann ${ref}.bwt ${ref}.pac ${ref}.sa

# Clean up
if ${cleanup}; then echo "Cleaning up input and intermediate files"
	rm ${ref}; rm ${ref}.fai; rm ${rb}.dict; rm ${ref}.amb; rm ${ref}.ann; rm ${ref}.bwt; rm ${ref}.pac; rm ${ref}.sa
else
	echo "Skipping clean up"
fi

