#!/usr/bin/env python
import re
import argparse
__author__ = 'Sama'
 
parser = argparse.ArgumentParser(description='This is a script to convert the alleles in the blocks of the samtools phased output into fasta format.\nThis fasta files are not separated in cymt numt yet.\nIt will give 2 output files; [output]_phase0.txt and [output]_phase1.txt')
parser.add_argument('-i','--input', help='Input fasta file name',required=True)
parser.add_argument('-s','--sam', help='Input samtools phased output file name',required=True)
parser.add_argument('-o','--output',help='Output file name', required=True)
args = parser.parse_args()
 
DICT={}

# Save in memory positions of phased blocks and alleles

with open(args.sam,"r") as file1:
	for line in file1:
		line=line.rstrip("\n")
		if re.match("PS",line):
			linesplit=line.split("\t")
			if linesplit[2] != linesplit[3]:
				name=str(linesplit[1]) + ":" + str(linesplit[2])
				DICT[name]={}
		
		elif re.match("M1",line):
			linesplit=line.split("\t")
			DICT[name][str((int(linesplit[3])-int(linesplit[2])))]=[linesplit[4],linesplit[5]]

# Make the two phased alleles in fasta format

with open(args.input,"r") as file1:
	z=1
	with open(args.output + ".phase0.txt","w") as my_file1:
		with open(args.output + ".phase1.txt","w") as my_file2:
# For every fasta block ...
			for line in file1:
				line=line.rstrip("\n")
				if re.match(">",line):
					linesplit=line.split(">")
					linelast=linesplit[1].split(":")
				else:
					chunk=str(linelast[0]) + ":" + str(linelast[1])
					phase0=list(line)
					phase1=list(line)
# ... two phased sequences are going to be formed with the phased snps
					for key in DICT[chunk]:
						phase0[int(key)] = DICT[chunk][key][0]
						phase1[int(key)] = DICT[chunk][key][1]
					phase0="".join(phase0)
					phase1="".join(phase1)
# Save output
					my_file1.write(">" + args.output + "_phase0_" + str(z) + "\n")
					my_file1.write(phase0 + "\n")
					my_file2.write(">" + args.output + "_phase1_" + str(z) + "\n")
					my_file2.write(phase1 + "\n")
					z+=1


