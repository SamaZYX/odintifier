#!/usr/bin/env bash

############################
#ADD HOW YOU CALL GATK HERE#
############################

#example GenomeAnalysisTK="/home/usr/bin/GenomeAnalysisTK" or simply "GenomeAnalysisTK"
GenomeAnalysisTK=""

############################

path_odintifier=$(dirname "${BASH_SOURCE[0]}")

if [ "$GenomeAnalysisTK" == "" ] ; then
	echo "Add the path of GATK and odintifier in this script"
	exit
fi

usage()
{
cat << EOF
usage: $0 options

This script phase with GATK and returns two consensus sequences and a bed file.

OPTIONS:
	-h|--help		Show this message
	-r|--ref		Reference sequence used in the mapping
	-b|--bam		BAM file (preferably without duplicates and realigned)
	-v|--vcf		VCF file
	-c|--compref		Sequence to compare the phasing blocks
	-o|--outbase		Output prefix
EOF
}

PARSED_OPTIONS=$(getopt -n "$0" -o hr:b:v:c:o: --long "help,ref:,bam:,vcf:,comref:,outbase:"  -- "$@")
 
if [ $? -ne 0 ]; then exit 1; fi

eval set -- "$PARSED_OPTIONS"

ref=""
bam=""
vcf=""
compref=""
outbase=""

while true; do
	case $1 in
		-h|--help) usage; exit 1 ;;
		-r|--ref) shift; ref=$1 ;;
		-b|--bam) shift; bam=$1 ;;
		-v|--vcf) shift; vcf=$1 ;;
		-c|--compref) shift; compref=$1 ;;
		-o|--outbase) shift; outbase=$1 ;;
		--) shift; break ;;
		*) break ;;
	esac
	shift
done

if [[ -z $ref ]] || [[ -z $bam ]] || [[ -z $vcf ]] || [[ -z $compref ]] || [[ -z $outbase ]] ; then usage; exit 1 ; fi

cymt_numt_gatk="${path_odintifier}/cymt_numt_gatk.py"
phased_consensus="${path_odintifier}/phased_consensus.R"

if [[ ! -x $cymt_numt_gatk ]] || [[ ! -x $phased_consensus ]] ; then 
	echo "Programs were not found or are not executable"
	echo "Check $cymt_numt_gatk or $phased_consensus"
	exit
fi

$GenomeAnalysisTK -T ReadBackedPhasing -R $ref -I $bam --variant $vcf -o ${vcf%.*}.phase.vcf --phaseQualityThresh 10.0

$cymt_numt_gatk -i $ref -v ${vcf%.*}.phase.vcf -o $outbase

Rscript $phased_consensus $compref ${outbase}.phase0.txt ${outbase}.phase1.txt ${outbase}.bed $ref $outbase



