#!/bin/bash

#SBATCH -J test_savp
#SBATCH -o test_savp.o%j
#SBATCH -A iPlant-Collabs
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 24

#---------------------------------------------
# Script for variant calling 
#
# Original author: marlies dolezal (may-june 2013)
# Alterations: christine baes (june, july, august 2013)
# Modifications: eric fritz-waters (september, october 2013)
# Testing mods for tacc systems:james carson (january 2015 - march 2016)
#---------------------------------------------

module load python/2.7.11
module load gatk/3.5.0

#--------------------------------------------
# Inputs and parameters
#--------------------------------------------

# The bam file to analyze - required
# Warning - .bam.bai file must also be in same directory
FILE="SRR2601715.mem.bam"

# The reference files - required
# Warning - must be the same as used to create the bam
# Warning - bundle must be have common base name as contents
# Warning - contents must include .dict   .fa    .fa.fai    .fa.{amb, ann, bwt, pac, sa}
reference_bundle="e-coli-K-12.tar"  

# Known variant vcf file - optional
dbSNP_work="known.variants.vcf"
if [ "${dbSNP_work}" = "" ]; then 
	realign_dbsnp_arg=""
	ug_dbsnp_arg=""
	run_recalibrator=false
	recal=""
else
	realign_dbsnp_arg="-known ${dbSNP_work}"
	ug_dbsnp_arg="--dbsnp ${dbSNP_work}"
	run_recalibrator=true
	recal=".recal"
fi

# Region (aka location) - optional
location=""
if [ "${location}" = "" ]; then 
	gatk_region_arg=""
	mpileup_region_arg=""
	region=""
else 
	gatk_region_arg="-L ${location}"
	mpileup_region_arg="-r ${location}"
	region=".${location}"
fi

run_unifiedgenotyper=true
run_mpileup=true
run_platypus=true

cleanup_inputs=false
cleanup_intermediates=true

# number of cores to use per function
n=44

# memory per function
m=48

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

tar -xvf bin.tgz

JAVADIR=
PICARD=./bin/picard-tools-1.141/
GATK=$TACC_GATK_DIR
SAMTOOLS=./bin/samtools-0.1.19/
BCFTOOLS=./bin/samtools-0.1.19/bcftools/
vcfutXXX=./bin/samtools-0.1.19/bcftools/vcfutXXX.pl 
PLATYPUS=./bin/Platypus_0.5.2/

#---------------------------------------------
# FUNCTIONS
#---------------------------------------------

function GATKpipeline () {
	runTimeLog=${log}${bull}${region}.runtime.txt
	RIGHT_NOW=$(date +"%s")
	echo "Beginning: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
	
	function extractRegion () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME " $RIGHT_NOW >> ${runTimeLog}
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T PrintReads --unsafe -nct ${n}  --defaultBaseQualities 30 \
		-R ${reference} \
		-I ${bull}.bam \
		-o ${bull}${region}.bam \
		${gatk_region_arg} \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME " ${RIGHT_NOW} >> ${runTimeLog}
	}
	if [ "${region}" != "" ]; then extractRegion; fi

	function RealignerTargetCreator () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" $RIGHT_NOW >> ${runTimeLog}
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T RealignerTargetCreator --unsafe -nt ${n} \
		-R ${reference} \
		-I ${bull}${region}.bam \
		${realign_dbsnp_arg} \
		-o ${temp}${bull}${region}.RTCintervals.list \
		${gatk_region_arg} \
		-l INFO \
		-log ${log}${bull}${region}.RealignerTargetCreator.log \
		-allowPotentiallyMisencodedQuals
		#attention -o file must be named .list or .bed for IndelRealigner
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	RealignerTargetCreator

	function IndelRealigner () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T IndelRealigner --unsafe \
		-R ${reference} \
		-I ${bull}${region}.bam \
		-targetIntervals ${temp}${bull}${region}.RTCintervals.list \
		${realign_dbsnp_arg} \
		${gatk_region_arg} \
		--consensusDeterminationModel USE_SW \
		-o ${temp}${bull}${region}.realigned.bam \
		-l INFO \
		-log ${log}${bull}${region}.IndelRealigner.log \
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
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T BaseRecalibrator --unsafe -nct ${n} \
		-R ${reference} \
		${gatk_region_arg} \
		-I ${temp}${bull}${region}.realigned.bam \
	     	-knownSites ${dbSNP_work} \
		-o ${temp}${bull}${region}.realigned.grp \
		-l INFO \
		-allowPotentiallyMisencodedQuals 
		#--plot_pdf_file ${Reg}${bull}.pre_recal.pdf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 1" ${RIGHT_NOW} >> ${runTimeLog}
		
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 2" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T BaseRecalibrator --unsafe -nct ${n} \
		-R ${reference} \
		${gatk_region_arg} \
		-I ${temp}${bull}${region}.realigned.bam \
		-BQSR ${temp}${bull}${region}.realigned.grp \
     		-knownSites ${dbSNP_work} \
		-o ${temp}${bull}${region}.realigned${recal}.grp \
		-l INFO \
		-allowPotentiallyMisencodedQuals #\
		#--plot_pdf_file ${Reg}${bull}.post_recal.pdf
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 2" ${RIGHT_NOW} >> ${runTimeLog}
			
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME Part 3" ${RIGHT_NOW} >> ${runTimeLog}		
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T PrintReads --unsafe -nct ${n} \
		-R ${reference} \
		${gatk_region_arg} \
		-I ${temp}${bull}${region}.realigned.bam \
		-BQSR ${temp}${bull}${region}.realigned.grp \
		-o ${temp}${bull}${region}.realigned${recal}.bam \
		-l INFO \
		-log ${log}${bull}${region}.realigned${recal}.bam.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME Part 3" ${RIGHT_NOW} >> ${runTimeLog}
	}
	if $run_recalibrator; then BaseRecalibrator; fi 

	function UnifiedGenotyper () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
		${JAVADIR}java -Xmx${m}g -jar ${GATK}GenomeAnalysisTK.jar \
		-T UnifiedGenotyper --unsafe -nt ${n} \
		-R ${reference} \
		${gatk_region_arg} \
		-I ${temp}${bull}${region}.realigned${recal}.bam \
		-o ${results}${bull}${region}.realigned${recal}.UG.raw.vcf \
		${ug_dbsnp_arg} \
		-out_mode EMIT_VARIANTS_ONLY \
		-stand_call_conf 30.0 \
		-stand_emit_conf 30.0 \
		--genotype_likelihoods_model BOTH \
		-l INFO  \
		-log ${log}${bull}${region}.realigned${recal}.UnifiedGenotyper.log \
		-allowPotentiallyMisencodedQuals
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME" ${RIGHT_NOW} >> ${runTimeLog}
	}
	if ${run_unifiedgenotyper}; then UnifiedGenotyper; fi

	function platypus () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
		python ${PLATYPUS}Platypus.py callVariants --nCPU=${n} --refFile=${reference} --bamFiles=${temp}${bull}${region}.realigned${recal}.bam --output=${results}${bull}${region}.realigned${recal}.platypus.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	if ${run_platypus}; then platypus; fi	
	
	function mpileup () {
		RIGHT_NOW=$(date +"%s")
		echo "BEGINNING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
		echo "${SAMTOOLS}samtools mpileup -f ${reference} ${mpileup_region_arg} -A -C50 -E -D -S -u ${temp}${bull}${region}.realigned${recal}.bam"
		${SAMTOOLS}samtools mpileup -f ${reference} ${mpileup_region_arg} -A -C50 -E -D -S -u ${temp}${bull}${region}.realigned${recal}.bam | \
		${BCFTOOLS}bcftools view -cvgb - > ${temp}${bull}${region}.raw3.bcf 2> ${temp}${bull}${region}.raw3.bcf.log 
		${BCFTOOLS}bcftools view ${temp}${bull}${region}.raw3.bcf | ${vcfutXXX} varFilter -D800 - > ${results}${bull}${region}.realigned${recal}.flt.vcf 
		RIGHT_NOW=$(date +"%s")
		echo "ENDING: $FUNCNAME 3rd" ${RIGHT_NOW} >> ${runTimeLog}
	}
	if ${run_mpileup}; then mpileup; fi
	
	RIGHT_NOW=$(date +"%s")
	echo "ENDING: $FUNCNAME " ${RIGHT_NOW} >> ${runTimeLog}
}
GATKpipeline

# Clean up
if ${cleanup_inputs}; then echo "Cleaning up input files"
	rm ${FILE}; rm ${FILE}.bai										# remove inputs
	rm ${dbSNP_work}; rm ${dbSNP_work}.idx          # remove dbSNP
	rm ${reference_bundle}; 
	rm ${reference}; rm ${dictionary} 		# remove reference files
	rm ${reference}.fai; rm ${reference}.amb; rm ${reference}.ann; 	# remove reference files
	rm ${reference}.bwt; rm ${reference}.pac; rm ${reference}.sa	# remove reference files
else
	echo "Skipping input file clean up"
fi

if ${cleanup_intermediates}; then echo "Cleaning up intermediate files"
	rm -rf ${temp}			# remove intermediate files
else
	echo "Skipping intermediate file clean up"
fi

rm -rf bin
#rm -rf bin.tgz
