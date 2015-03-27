#!/usr/bin/env bash

path_odintifier=$(dirname "${BASH_SOURCE[0]}")

usage()
{
cat << EOF
usage: $0 options

This script phase with GATK and returns two consensus sequences and a bed file.

OPTIONS:
	-h|--help		Show this message
	-r|--ref		Reference sequence used in the mapping
	-0|--bam0		BAM phase 0 from samtools
	-1|--bam1		BAM phase 1 from samtools
	-p|--phase		Phasing output from samtools
	-b|--blocks		Phasing blocks from the program phased_consensus_samtools.R
	-o|--outbase		Output prefix
EOF
}

PARSED_OPTIONS=$(getopt -n "$0" -o hr:0:1:p:b:o: --long "help,ref:,bam0:,bam1:,phase:,blocks:,outbase:"  -- "$@")
 
if [ $? -ne 0 ]; then exit 1; fi

eval set -- "$PARSED_OPTIONS"

ref=""
bam0=""
bam1=""
phase=""
blocks=""
outbase=""

while true; do
	case $1 in
		-h|--help) usage; exit 1 ;;
		-r|--ref) shift; ref=$1 ;;
		-0|--bam0) shift; bam0=$1 ;;
		-1|--bam1) shift; bam1=$1 ;;
		-p|--phase) shift; phase=$1 ;;
		-b|--blocks) shift; blocks=$1 ;;
		-o|--outbase) shift; outbase=$1 ;;
		--) shift; break ;;
		*) break ;;
	esac
	shift
done

if [[ -z $ref ]] || [[ -z $bam0 ]] || [[ -z $bam1 ]] || [[ -z $phase ]] || [[ -z $blocks ]] || [[ -z $outbase ]] ; then usage; exit 1 ; fi

fixing_samtools_output_R="${path_odintifier}/fixing_samtools_output.R"

if [[ ! -x $fixing_samtools_output_R ]] ; then 
	echo "Program was not found or is not executable"
	exit
fi


# Start fixing

awk '{if($0~/^M1/){x=$4-1;print $2"\t"x"\t"$4}}' $phase > ${outbase}.phase.M1.bed
bedtools getfasta -fi $ref -bed ${outbase}.phase.M1.bed -fo ${outbase}.ref.M1.fasta

#bam0
samtools mpileup -uf $ref $bam0 | bcftools view -cg - | vcfutils.pl vcf2fq | awk 'BEGIN{flag=0}{if(flag==0){if($0~/^+/){flag=1}else{if($0~/^@/){sub(/@/,">",$0);print $0}else{print $0}}}}' > ${outbase}.phase.bam0.fasta
bedtools getfasta -fi ${outbase}.phase.bam0.fasta -bed ${outbase}.phase.M1.bed -fo ${outbase}.phase.bam0.M1.fasta
#bam1
samtools mpileup -uf $ref $bam1 | bcftools view -cg - | vcfutils.pl vcf2fq | awk 'BEGIN{flag=0}{if(flag==0){if($0~/^+/){flag=1}else{if($0~/^@/){sub(/@/,">",$0);print $0}else{print $0}}}}' > ${outbase}.phase.bam1.fasta
bedtools getfasta -fi ${outbase}.phase.bam1.fasta -bed ${outbase}.phase.M1.bed -fo ${outbase}.phase.bam1.M1.fasta

egrep -e "^(PS|M1)" $phase > ${outbase}.PS_M1.phase.out

Rscript $fixing_samtools_output_R $blocks ${outbase}.PS_M1.phase.out ${outbase}.phase.bam0.M1.fasta ${outbase}.phase.bam1.M1.fasta ${outbase}

