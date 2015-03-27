#!/usr/bin/python
import re
import argparse
__author__ = 'Sama'
 
parser = argparse.ArgumentParser(description='This is a demo script to retrieve the cymt and numt using samtools output by Sama.')
parser.add_argument('-i','--input', help='Input fasta file name',required=True)
parser.add_argument('-s','--sam', help='Input samtools phase output file name',required=True)
parser.add_argument('-o','--output',help='Output file name', required=True)
args = parser.parse_args()
 
DICT={}

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


with open(args.input,"r") as file1:
	z=1
	with open(args.output + ".phase0.txt","w") as my_file1:
		with open(args.output + ".phase1.txt","w") as my_file2:
			for line in file1:
				line=line.rstrip("\n")
				if re.match(">",line):
					linesplit=line.split(">")
					linelast=linesplit[1].split(":")
				else:
					chunk=str(linelast[0]) + ":" + str(linelast[1])
					father=list(line)
					mother=list(line)
					for key in DICT[chunk]:
						father[int(key)] = DICT[chunk][key][0]
						mother[int(key)] = DICT[chunk][key][1]
					father="".join(father)
					mother="".join(mother)
					my_file1.write(">" + args.output + "_phase0_" + str(z) + "\n")
					my_file1.write(father + "\n")
					my_file2.write(">" + args.output + "_phase1_" + str(z) + "\n")
					my_file2.write(mother + "\n")
					z+=1


