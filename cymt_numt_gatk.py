#!/usr/bin/env python
import re
import argparse
__author__='Sama'

parser=argparse.ArgumentParser(description='This is a script to convert the alleles in the blocks of the phased vcf into fasta format.\nThis fasta files are not separated in cymt numt yet.\nIt will give 3 output files; [output]_phase0.txt, [output]_phase1.txt and [output].bed')
parser.add_argument('-i','--input',help='Input fasta file name',required=True)
parser.add_argument('-v','--vcf',help='Input phased vcf file name',required=True)
parser.add_argument('-o','--output',help='Output file name',required=True)
args=parser.parse_args()

# Save referene fasta in memory

refseq=""
with open(args.input,"r") as file1:
	for line in file1:
		line=line.rstrip("\n")
		if not re.match(">",line):
			refseq+=line
refseq=list(refseq)

bed=open("%s.bed"%(args.output),"w")

seqphase0=[]
seqphase1=[]
seqphase0.append("")
seqphase1.append("")
i=0
j=0
nv=0
nv2=0
start=1
with open(args.vcf,"r") as file1:
	for line in file1:
		line=line.rstrip("\n")
		if not re.match("#",line):
			linesplit=line.split("\t")
			linelast=linesplit[9].split(":")
# Only analize biallelic snps
			if len(linesplit[3]) == 1 and len(linesplit[4])==1:
# Make the proper allele of every block depending on the phase tag in the vcf 
				if int(linesplit[1])-1 != j:
					seqphase1[i]+="".join(list(refseq[j:int(linesplit[1])-1]))
					seqphase0[i]+="".join(list(refseq[j:int(linesplit[1])-1]))
					j=int(linesplit[1])-1
				j+=1
# Starting a phased block or single snp
				if"0/1" == linelast[0]:
					if"LowQual" == linesplit[6]:
						seqphase0[i]+=linesplit[3]
						seqphase1[i]+=linesplit[3]
					else:
						if nv2 < 2:
							seqphase0[i]=""
							seqphase1[i]=""
						else:
							if int(linesplit[1])-1 != nv:
								seqphase0[i]=seqphase0[i][:-(int(linesplit[1])-nv-1)]
								seqphase1[i]=seqphase1[i][:-(int(linesplit[1])-nv-1)]
							seqphase0.append("")
							seqphase1.append("")
							i+=1
							bed.write("%s\t%s\t%s\n" %(linesplit[0],start-1,nv))
						start=int(linesplit[1])
						nv=int(linesplit[1])
						seqphase0[i]+=linesplit[3]
						seqphase1[i]+=linesplit[4]
						nv2=1
# Starting a phased block or homozygous snp
				elif"1/1"==linelast[0]:
					if"LowQual" == linesplit[6]:
						seqphase0[i]+=linesplit[3]
						seqphase1[i]+=linesplit[3]
					else:
						if nv2 < 2:
							seqphase0[i]=""
							seqphase1[i]=""
						else:
							if int(linesplit[1])-1 != nv:
								seqphase0[i]=seqphase0[i][:-(int(linesplit[1])-nv-1)]
								seqphase1[i]=seqphase1[i][:-(int(linesplit[1])-nv-1)]
							seqphase0.append("")
							seqphase1.append("")
							i+=1
							bed.write("%s\t%s\t%s\n" %(linesplit[0],start-1,nv))
						start=int(linesplit[1])
						nv=int(linesplit[1])
						seqphase0[i]+=linesplit[4]
						seqphase1[i]+=linesplit[4]
						nv2=1
# Phased heterozygous snp
				elif"1|0" == linelast[0]:
					seqphase0[i]+=linesplit[4]
					seqphase1[i]+=linesplit[3]
					nv=int(linesplit[1])
					nv2+=1
# Phased heterozygous snp
				elif"0|1" == linelast[0]:
					seqphase0[i]+=linesplit[3]
					seqphase1[i]+=linesplit[4]
					nv=int(linesplit[1])
					nv2+=1
# Phased homozygous snp
				elif"1|1" == linelast[0]:
					seqphase0[i]+=linesplit[4]
					seqphase1[i]+=linesplit[4]
					nv=int(linesplit[1])
					nv2+=1
				else:
					print"weirdstuff"+linelast[0]+"\n"+line+"\n"
#For the last snp in the vcf file	
if nv2 < 2:
	seqphase0.pop()
	seqphase1.pop()
else:
	bed.write("%s\t%s\t%s\n" %(linesplit[0],start-1,nv))

bed.close()

# Write the phased alleles in two files

z=1
with open(args.output + ".phase0.txt","w") as my_file:
	for j in seqphase0:
		my_file.write(">" + args.output + "_phase0_" + str(z) + "\n")
		z+=1
		my_file.write(j + "\n") 

z=1
with open(args.output + ".phase1.txt","w") as my_file:
	for j in seqphase1:
		my_file.write(">" + args.output + "_phase1_" + str(z) + "\n")
		z+=1
		my_file.write(j + "\n")

