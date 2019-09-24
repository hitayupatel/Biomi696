#! /usr/bin/env python

import Bio.Seq
from Bio import SeqIO


def file_read(filename):
    #with open(filename, 'r') as infile:
	read_length=0
	lines_count=0
	for record in SeqIO.parse(filename, "fastq"):    
		read_length += len(record)
		lines_count += 1
	return read_length/lines_count
	

def main(file1, file2):
    file1_average = file_read(file1)
    file2_average = file_read(file2)
	print((file1_average + file2_average)/2)


if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2])