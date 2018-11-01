#!/usr/bin/env python

import argparse
import gzip

#*******************************************************************************
#Function	: get_arguments
#Description: get the arguments from the command line for generalization 
#				of the program
#Parameters	: none
#Returned	: parse_args - the arguments needed for the program
#******************************************************************************* 

def get_arguments():
	parser = argparse.ArgumentParser (description = "dedupe reads in a sam file")
	parser.add_argument("-f", "--file", help="full absolute path to sam file to dedupe ", required=True, type=str)
	parser.add_argument("-p", "--paired", help="is the file paried end??", required=False, action='store_true')
	parser.add_argument("-u", "--umi", help="file containing UMIs used in your reads", required=False, type=str)
	return parser.parse_args()
args = get_arguments()

#*******************************************************************************
#Function	  : writeToFile
#Description  : write unique reads to a file
#Parameters	  : read - unique read to be written to a file
#               file - output file
#Returned	  : none
#*******************************************************************************

def writeToFile(out, dictionary):
	for value in dictionary.values():
		out.write(str(value))
		out.write("\n")
		
#*******************************************************************************
#Function	  : adjustPosForward
#Description  : adjust the left most position taking into account soft clipping 
#Parameters	  : read - sam file read with information on soft clipping and position
#Returned	  : position - adjusted position of read based on CIGAR string
#*******************************************************************************
import re

def adjustPosForward(read):
	adjustment = 0
	unit = ""
	first = True
	read = read.strip()
	read = read.split("\t")
	position = int(read[3])
	cigar = read[5]
	for char in cigar:
		if str.isdigit(char): 
			unit = unit + char
		else:
			unit = unit + char
			#adjust for left S
			mMatch = re.match("([0-9]+)S", unit)
			if mMatch:
				adjust = mMatch.group(1)
				if first == True:
					adjustment = adjustment - int(adjust)
			break
	return (position + adjustment)

#*******************************************************************************
#Function	  : adjustPosReverse
#Description  : adjust the left most position to create right most mapping position 
#              taking into account soft clipping 
#Parameters	  : read - sam file read with information on soft clipping and position
#Returned	  : position - adjusted position of read based on CIGAR string
#*******************************************************************************

def adjustPosReverse(read):
	adjustment = 0
	unit = ""
	first = True
	read = read.strip()
	read = read.split("\t")
	position = int(read[3])
	cigar = read[5]
	for char in cigar:
		if str.isdigit(char): 
			unit = unit + char
		else:
			unit = unit + char
			#adjust for left S IGNORE
			#if first:
				#adjustment = adjustment #+ adjustPos("S", unit)
			#adjust for right S
			if not first:
				adjustment = adjustment + adjustPos("S", unit)
			#adjust for M
			adjustment = adjustment + adjustPos("M", unit)
			#adjust for I IGNORE
			#adjustment = adjustment #- adjustPos("I", unit)
			#adjust for D
			adjustment = adjustment + adjustPos("D", unit)
			#adjust for N
			adjustment = adjustment + adjustPos("N", unit)
			unit = ""
			first = False
	return position + (adjustment - 1)

#*******************************************************************************
#Function	  : adjustPos
#Description  : adjust the left most position to create right most mapping position 
#              taking into account soft clipping 
#Parameters	  : letter - what condtion is being tested for in the cigar string
#             : unit   - one number and letter unit from the cigar string
#Returned	  : int    - number that the position needs to be changed by
#*******************************************************************************

def adjustPos(letter, unit):
	condition = re.match("([0-9]+)%s" % letter, unit)
	if condition:
		return int(condition.group(1))
	else:
		return 0
	
#*******************************************************************************
#Function	  : getUMI
#Description  : get the UMI from the header line (last 8 characters in col 0) 
#Parameters	  : read - sam file read with header column
#Returned	  : UMI - the UMI for the read
#*******************************************************************************

def getUMI(read):
	read = read.strip()
	read = read.split("\t")
	read = read[0]
	UMI = read[-8:]
	return UMI

#*******************************************************************************
#Function	  : getChrom
#Description  : get the chromosome from the read (col 2) 
#Parameters	  : read - sam file read with header column
#Returned	  : chrom - the chromosome for the read
#*******************************************************************************

def getChrom(read):
	read = read.strip()
	read = read.split("\t")
	return read[2]

#*******************************************************************************
#Function	  : reverseStrand
#Description  : get the strandedness from bitwise flag (col 1) 
#Parameters	  : read - sam file read 
#Returned	  : bool - true if reverse strand
#*******************************************************************************

def reverseStrand(read):
	read = read.strip()
	read = read.split("\t")
	read = int(read[1])
	return(read & 16 == 16)

#from subprocess import call

#for command in ("samtools view -bS args.file | samtools sort -O sam -args.file_sorted"):
	# shell=True is so you can handle redirects like in the 3rd command
#    call(command, shell=True)
	
if args.umi:
	UMI_list = True
	UMIs = {}
	with open(args.umi, "r") as fh:
		for line in fh:
			line = line.strip()
			UMIs[line] = ""
else:
	UMI_list = False
	print("Continuing with no UMI reference")
	
			
if args.paired:
	pairedEnd = True
else:
	pairedEnd = False
	print("Continuing as single end reads")
	
			
counter = 0
chromosome = 1
unique = {}

with open(args.file, "r") as fh:
	with open(args.file + "_deduped", "w") as out:
		for line in fh:
			counter = counter + 1
			if counter % 100000 == 0:
				print("On line:", counter)
			if line.startswith("@"):
				out.write(line)
			else:
				line = line.strip()
				if chromosome != getChrom(line):
					chromosome = getChrom(line)
					writeToFile(out, unique)
					unique = {}
				UMI = getUMI(line)
				if (UMI_list == True and UMI in UMIs) or UMI_list == False:
					if pairedEnd == True:
						if reverseStrand(line) == True:
							position = adjustPosForward(line)
						else:
							position = adjustPosForward(line)
					else:
						position = adjustPosForward(line)
					identifier = UMI + str(position)
					if identifier not in unique:
						unique[UMI, position] = line
		writeToFile(out, unique)