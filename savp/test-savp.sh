#!/bin/bash

#SBATCH -J test_savp
#SBATCH -o savp_test.o%j
#SBATCH -A iPlant-Collabs
#SBATCH -p normal
#SBATCH -t 47:30:00
#SBATCH -N 1
#SBATCH -n 16

#---------------------------------------------
# Script for variant calling 
#
# Original author: marlies dolezal (may-june 2013)
# Alterations: christine baes (june, july, august 2013)
# Modifications: eric fritz-waters (september, october 2013)
# Testing mods for stampede:	james carson (january 2015)
#---------------------------------------------

#--------------------------------------------
# Inputs and parameters
#--------------------------------------------

# The bam file to analyze
# Warning - .bam.bai file must also be in same directory
FILE="JER.33680.AP.01.mem.fixed.bam"

# The reference files.  
# Warning - must be the same as used to create the bam
# Warning - bundle must be have common base name as contents
# Warning - contents must include .dict   .fa    .fa.fai    .fa.{amb, ann, bwt, pac, sa}
reference_bundle="umd_3_1_Y_Mito.tar"  

# Known variant vcf file
dbSNP_work="Bos_taurus.dbSNP.vcf"

#---------------------------------------------
# DIRECTORIES
#---------------------------------------------

bull=`echo "$FILE" | sed 's/.bam//g'`
temp=${bull}.temp/
log=${bull}.logs/
results=${bull}.results/

mkdir $temp
mkdir $log
mkdir $results

#---------------------------------------------
# PARAMETER FILES
#---------------------------------------------

# eventually update the following lines to access other types of bundles
tar -xvf ${reference_bundle}
ref=`echo "$reference_bundle" | sed 's/.tar//g'`
reference=${ref}.fa
dictionary=${ref}.dict


#---------------------------------------------
# PROGRAMME VERSIONS
#---------------------------------------------

# tar -xvf bin.tgz
# export PATH=$PATH:"$PWD/bin"
#
# then at end add
# rm -rf bin



JAVADIR=

module load picard
PICARD=$TACC_PICARD_DIR

module load gatk/3.2.2
GATK=$TACC_GATK_DIR

module load samtools
SAMTOOLS=
BCFTOOLS=
vcfutils=vcfutils.pl # comes with samtools module

module load python
module load Platypus
platypus=$TACC_PLATYPUS_DIR



#---------------------------------------------
# FUNCTIONS
#---------------------------------------------

function GATKpipeline () {
	runTimeLog=${log}${bull}.runtime.txt
	RIGHT_NOW=$(date +"%s")
	echo "Beginning: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
	
	function RealignerTargetCreator () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
		${JAVADIR}java -Xmx26g -jar ${GATK}GenomeAnalysisTK.jar \
		-T RealignerTargetCreator --unsafe -nt 16 \
		-R ${reference} \
		-I ${bull}.bam \
		-known ${dbSNP_work} \
		-o ${temp}${bull}.RTCintervals.list \
		-l INFO \
		-log ${log}${bull}.RealignerTargetCreator.log \
		-allowPotentiallyMisencodedQuals
		#attention -o file must be named .list or .bed for IndelRealigner
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	RealignerTargetCreator

	function IndelRealigner () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx26g -jar ${GATK}GenomeAnalysisTK.jar \
		-T IndelRealigner --unsafe \
		-R ${reference} \
		-I ${bull}.bam \
		-targetIntervals ${temp}${bull}.RTCintervals.list \
		-known ${dbSNP_work} \
		--consensusDeterminationModel USE_SW \
		-o ${temp}${bull}.realigned.bam \
		-l INFO \
		-log ${log}${bull}.IndelRealigner.log \
		-allowPotentiallyMisencodedQuals
		#@ -known Non-indel variants in these files will be ignored. so we can use
		#batch download of all variants from NCBI or ensembl
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	IndelRealigner

	function BaseRecalibrator () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 1" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx26g -jar ${GATK}GenomeAnalysisTK.jar \
		-T BaseRecalibrator --unsafe -nct 16 \
		-R ${reference} \
		-I ${temp}${bull}.realigned.bam \
	     	-knownSites ${dbSNP_work} \
		-o ${temp}${bull}.realigned.grp \
		-l INFO \
		-allowPotentiallyMisencodedQuals 
		#--plot_pdf_file ${Reg}${bull}.pre_recal.pdf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 1" ${RIGHT_NOW} >> ${runTimeLog}
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 2" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx26g -jar ${GATK}GenomeAnalysisTK.jar \
		-T BaseRecalibrator --unsafe -nct 16 \
		-R ${reference} \
		-I ${temp}${bull}.realigned.bam \
		-BQSR ${temp}${bull}.realigned.grp \
     		-knownSites ${dbSNP_work} \
		-o ${temp}${bull}.realigned.recal.grp \
		-l INFO \
		-allowPotentiallyMisencodedQuals #\
		#--plot_pdf_file ${Reg}${bull}.post_recal.pdf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 2" ${RIGHT_NOW} >> ${runTimeLog}
			
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 3" ${RIGHT_NOW} >> ${runTimeLog}		
		${JAVADIR}java -Xmx26g -jar ${GATK}GenomeAnalysisTK.jar \
		-T PrintReads --unsafe -nct 16 \
		-R ${reference} \
		-I ${temp}${bull}.realigned.bam \
		-BQSR ${temp}${bull}.realigned.grp \
		-o  ${temp}${bull}.realigned.recal.bam \
		-l INFO \
		-log ${log}${bull}.realigned.recal.bam.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 3" ${RIGHT_NOW} >> ${runTimeLog}
	}
	BaseRecalibrator

	function UnifiedGenotyper () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx26g -jar ${GATK}GenomeAnalysisTK.jar \
		-T UnifiedGenotyper --unsafe -nt 16 \
		-R ${reference} \
		-I ${temp}${bull}.realigned.recal.bam \
		-o ${results}${bull}.UG.raw.vcf \
		--dbsnp $dbSNP_work \
		-out_mode EMIT_VARIANTS_ONLY \
		-stand_call_conf 30.0 \
		-stand_emit_conf 30.0 \
		--genotype_likelihoods_model BOTH \
		-l INFO  \
		-log ${log}${bull}.UnifiedGenotyper.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	UnifiedGenotyper

	function mpileup () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
		${SAMTOOLS}samtools mpileup -f ${reference} -A -C50 -E -D -S -u ${temp}${bull}.realigned.recal.bam | \
		${BCFTOOLS}bcftools view -cvgb - > ${temp}${bull}.raw3.bcf 2> ${temp}${bull}.raw3.bcf.log 
		${BCFTOOLS}bcftools view ${temp}${bull}.raw3.bcf | ${vcfutils} varFilter -D800 - > ${results}${bull}.realigned.recal.flt.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	mpileup
		
	function platypus () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
		python ${platypus}Platypus.py callVariants --nCPU=16 --refFile=${reference} --bamFiles=${temp}${bull}.realigned.recal.bam --output=${results}${bull}.realigned.recal.platypus.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	platypus
	
	RIGHT_NOW=$(date +"%s")
	echo "ENDING: $FUNCNAME " ${RIGHT_NOW} >> ${runTimeLog}
}
GATKpipeline