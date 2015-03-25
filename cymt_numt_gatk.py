#!/usr/bin/env python
import re
import argparse
__author__='Sama'

parser=argparse.ArgumentParser(description='This is a demo script to retrieve the cymt and numt by Sama.\nIt will give 3 output files; [output]_father.txt, [output]_mother.txt and [output].bed')
parser.add_argument('-i','--input',help='Input fasta file name',required=True)
parser.add_argument('-v','--vcf',help='Input vcf file name',required=True)
parser.add_argument('-o','--output',help='Output file name',required=True)
args=parser.parse_args()

refseq=""
with open(args.input,"r") as file1:
	for line in file1:
		line=line.rstrip("\n")
		if not re.match(">",line):
			refseq+=line
refseq=list(refseq)

bed=open("%s.bed"%(args.output),"w")
					
seqfather=[]
seqmother=[]
seqfather.append("")
seqmother.append("")
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
			if len(linesplit[3]) == 1 and len(linesplit[4])==1:
				if int(linesplit[1])-1 != j:
					seqmother[i]+="".join(list(refseq[j:int(linesplit[1])-1]))
					seqfather[i]+="".join(list(refseq[j:int(linesplit[1])-1]))
					j=int(linesplit[1])-1
				j+=1
				if"0/1" == linelast[0]:
					if"LowQual" == linesplit[6]:
						seqfather[i]+=linesplit[3]
						seqmother[i]+=linesplit[3]
					else:
						if nv2 < 2:
							seqfather[i]=""
							seqmother[i]=""
						else:
							if int(linesplit[1])-1 != nv:
								seqfather[i]=seqfather[i][:-(int(linesplit[1])-nv-1)]
								seqmother[i]=seqmother[i][:-(int(linesplit[1])-nv-1)]
							seqfather.append("")
							seqmother.append("")
							i+=1
							bed.write("%s\t%s\t%s\n" %(linesplit[0],start-1,nv))
						start=int(linesplit[1])
						nv=int(linesplit[1])
						seqfather[i]+=linesplit[3]
						seqmother[i]+=linesplit[4]
						nv2=1
				elif"1/1"==linelast[0]:
					if"LowQual" == linesplit[6]:
						seqfather[i]+=linesplit[3]
						seqmother[i]+=linesplit[3]
					else:
						if nv2 < 2:
							seqfather[i]=""
							seqmother[i]=""
						else:
							if int(linesplit[1])-1 != nv:
								seqfather[i]=seqfather[i][:-(int(linesplit[1])-nv-1)]
								seqmother[i]=seqmother[i][:-(int(linesplit[1])-nv-1)]
							seqfather.append("")
							seqmother.append("")
							i+=1
							bed.write("%s\t%s\t%s\n" %(linesplit[0],start-1,nv))
						start=int(linesplit[1])
						nv=int(linesplit[1])
						seqfather[i]+=linesplit[4]
						seqmother[i]+=linesplit[4]
						nv2=1
				elif"1|0" == linelast[0]:
					seqfather[i]+=linesplit[4]
					seqmother[i]+=linesplit[3]
					nv=int(linesplit[1])
					nv2+=1
				elif"0|1" == linelast[0]:
					seqfather[i]+=linesplit[3]
					seqmother[i]+=linesplit[4]
					nv=int(linesplit[1])
					nv2+=1
				elif"1|1" == linelast[0]:
					seqfather[i]+=linesplit[4]
					seqmother[i]+=linesplit[4]
					nv=int(linesplit[1])
					nv2+=1
				else:
					print"weirdstuff"+linelast[0]+"\n"+line+"\n"
	
if nv2 < 2:
	seqfather.pop()
	seqmother.pop()
else:
	bed.write("%s\t%s\t%s\n" %(linesplit[0],start-1,nv))

bed.close()

z=1
with open(args.output + ".phase0.txt","w") as my_file:
	for j in seqfather:
		my_file.write(">" + args.output + "_phase0_" + str(z) + "\n")
		z+=1
		my_file.write(j + "\n")

z=1
with open(args.output + ".phase1.txt","w") as my_file:
	for j in seqmother:
		my_file.write(">" + args.output + "_phase1_" + str(z) + "\n")
		z+=1
		my_file.write(j + "\n")

