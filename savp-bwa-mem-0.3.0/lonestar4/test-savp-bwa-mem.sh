#!/bin/bash
#$ -V
#$ -cwd # Start job in submission directory
#$ -N test-savp_bwa_aln # Job Name
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12
#$ -q development       # Queue name normal
#$ -A iPlant-Collabs
#$ -l h_rt=01:00:00     # Run time (hh:mm:ss) 
#$ -M jcarson@tacc.utexas.edu   # Address for email notification
#$ -m be        # Email at Begin and End of job

module load jdk64
module load picard
module load samtools
module load bwa/0.7.7

## INPUTS
data_a="HO.2081.AP.03.1_Unique.fastq.gz"
data_b="HO.2081.AP.03.2_Unique.fastq.gz"
barcode="HO.2081.AP.03"   # make the user input this, required string
reference_bundle="umd_3_1_Y_Mito.tar"  
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

# number of cores to use per function
n=11

## MAIN FUNCTIONS


## If the reads are not paired, simply run samse, otherwise build the sam index for the mates and align them with sampe
if [ -z ${data_b} ]
    then bwa mem -t ${n} ${reference} ${data_a} | samtools view -bhS -@${n} - | samtools sort -@${n} - out_fasta_sorted
else
    bwa mem -t ${n} ${reference} ${data_a} ${data_b} | samtools view -bhS -@${n} - | samtools sort -@${n} - out_fasta_sorted
fi

java -Xmx12g -jar $TACC_PICARD_DIR/AddOrReplaceReadGroups.jar INPUT=out_fasta_sorted.bam OUTPUT=${barcode}.mem.bam RGLB=1 RGPL=illumina RGPU=all RGSM=${barcode} VALIDATION_STRINGENCY=SILENT

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
