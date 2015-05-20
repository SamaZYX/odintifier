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
	-c|--compref		Sequence to compare the phasing blocks
	-o|--outbase		Output prefix
	-t|--tmp		Leave the temporal files, don't delete them
EOF
}

PARSED_OPTIONS=$(getopt -n "$0" -o htr:b:c:o: --long "help,tmp,ref:,bam:,compref:,outbase:"  -- "$@")

if [ $? -ne 0 ]; then exit 1; fi

eval set -- "$PARSED_OPTIONS"

ref=""
bam=""
vcf=""
compref=""
outbase=""
tmp="1"

while true; do
	case $1 in
		-h|--help) usage; exit 1 ;;
		-t|--tmp) tmp="0";;
		-r|--ref) shift; ref=$1 ;;
		-b|--bam) shift; bam=$1 ;;
		-c|--compref) shift; compref=$1 ;;
		-o|--outbase) shift; outbase=$1 ;;
		--) shift; break ;;
		*) break ;;
	esac
	shift
done

if [[ -z $ref ]] || [[ -z $bam ]] || [[ -z $compref ]] || [[ -z $outbase ]] ; then usage; exit 1 ; fi

cymt_numt_samtools="${path_odintifier}/cymt_numt_samtools.py"
phased_consensus_samtools="${path_odintifier}/phased_consensus_samtools.R"
fixing_samtools_output="${path_odintifier}/fixing_samtools_output.sh"
fixing_samtools_output_R="${path_odintifier}/fixing_samtools_output.R"
phase_gatk="${path_odintifier}/phase_gatk.sh"

if [[ ! -x $cymt_numt_samtools ]] || [[ ! -x $phased_consensus_samtools ]] || [[ ! -x $fixing_samtools_output ]] || [[ ! -x $phase_gatk ]] || [[ ! -x $fixing_samtools_output_R ]] ; then 
	echo "Programs were not found or are not executable"
	echo -e "Check \n $cymt_numt_samtools \n $phased_consensus_samtools \n $fixing_samtools_output \n $phase_gatk \n $fixing_samtools_output_R"
	exit
fi


# Phase with samtools

bedtools genomecov -ibam $bam -bga > ${outbase}.DOC.txt
awk '{if($4<1){print $1"\t"$2"\t"$3"\t"$4 }}' ${outbase}.DOC.txt > ${outbase}.md1.bed

samtools calmd -AEur $bam $ref | samtools phase -Fb ${outbase}.phase - > ${outbase}.phase.out #Not fixing chimeric reads (F).
awk '{if($0~/^PS/){if(($4-$3)>0){print $2"\t"$3-1"\t"$4"\t"$2":"$3":"$4}}}' ${outbase}.phase.out > ${outbase}.phase.PS.bed
bedtools getfasta -fi $ref -bed ${outbase}.phase.PS.bed -name -fo stdout | awk '{if($0~/^>/){head = $0}else{if(length($0)>1){print head"\n"$0}}}' > ${outbase}.phase.PS.fasta

$cymt_numt_samtools -i ${outbase}.phase.PS.fasta -s ${outbase}.phase.out -o $outbase

Rscript $phased_consensus_samtools $compref ${outbase}.phase0.txt ${outbase}.phase1.txt ${outbase}.phase.PS.bed $ref $outbase

$fixing_samtools_output -r $ref -0 ${outbase}.phase.0.bam -1 ${outbase}.phase.1.bam -p ${outbase}.phase.out -b ${outbase}.choose_blocks.txt -o $outbase


# Making a new bam file

bedtools complement -i ${outbase}.phase.PS.bed -g ${ref}.fai > ${outbase}.complement.bed
bedtools intersect -abam ${outbase}.phase.0.bam -b ${outbase}.complement.bed -f 1.00 > ${outbase}.intersect.phase_0.bam
bedtools intersect -abam ${outbase}.phase.1.bam -b ${outbase}.complement.bed -f 1.00 > ${outbase}.intersect.phase_1.bam
samtools merge ${outbase}.intersect.merge.bam ${outbase}.intersect.phase_0.bam ${outbase}.intersect.phase_1.bam
samtools sort ${outbase}.intersect.merge.bam ${outbase}.intersect.merge.sort

num=0
while read f1
do
num=$[num+1]
head -n $num ${outbase}.phase.PS.bed | tail -n 1 > ${outbase}.block${num}.bed
bedtools intersect -abam ${outbase}.phase.${f1}.bam -b ${outbase}.block${num}.bed > good_block.${num}.${outbase}.intersect.bam
rm ${outbase}.block${num}.bed
done < ${outbase}.new_choose_blocks.txt

samtools merge ${outbase}.merge.all.bam ${outbase}.intersect.merge.sort.bam good_block.*.${outbase}.intersect.bam
samtools sort ${outbase}.merge.all.bam ${outbase}.merge.all.sort
samtools index ${outbase}.merge.all.sort.bam


# Phase with GATK and consensus

$GenomeAnalysisTK -T UnifiedGenotyper -R $ref -I ${outbase}.merge.all.sort.bam -stand_call_conf 30.0 -stand_emit_conf 10.0 -glm SNP -dcov 300 --out ${outbase}.samtools_gatk.vcf --output_mode EMIT_VARIANTS_ONLY

$phase_gatk -r $ref -b ${outbase}.merge.all.sort.bam -v ${outbase}.samtools_gatk.vcf -c $compref -o ${outbase}.phase_samtools


# Remove temporal files

if [ $tmp == "1" ] ; then
	rm ${outbase}.DOC.txt ${outbase}.md1.bed ${outbase}.phase.0.bam ${outbase}.phase.1.bam
	rm ${outbase}.phase.M1.bed ${outbase}.phase.PS.bed ${outbase}.phase.PS.fasta ${outbase}.phase0.txt
	rm ${outbase}.phase1.txt ${outbase}.choose_blocks.txt ${outbase}.complement.bed ${outbase}.intersect.phase_0.bam
	rm ${outbase}.intersect.phase_1.bam ${outbase}.intersect.merge.bam ${outbase}.intersect.merge.sort.bam good_block.*.${outbase}.intersect.bam
	rm ${outbase}.merge.all.bam ${outbase}.new_choose_blocks.txt ${outbase}.PS_M1.phase.out ${outbase}.ref.M1.fasta
	rm ${outbase}.phase.bam*.fasta*
fi

