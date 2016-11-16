#!/usr/bin/env python

"""Compare the results from multiple
PCR primers to determine which sets
would be the best for diagnostic assays"""

import glob
import sys
import os
import operator

results = {}

start_dir = os.getcwd()

def get_seq_name(fasta_in):
    name = os.path.basename(fasta_in)
    return name

results = {}

for infile in glob.glob(os.path.join(start_dir, '*_results.txt')):
    forward_dict = {}
    reverse_dict = {}
    name = get_seq_name(infile)
    reduced = name.replace("_results.txt", "")
    for line in open(infile, "U"):
        fields = line.split()
        info_fields= fields[0].split("_")
        direction = info_fields[1]
        if direction == "forward":
            try:
                forward_dict[reduced].append(fields[2])
            except KeyError:
                forward_dict[reduced]=[fields[2]]
        elif direction == "reverse":
            try:
                reverse_dict[reduced].append(fields[2])
            except KeyError:
                reverse_dict[reduced]=[fields[2]]
    for k,v in forward_dict.iteritems():
        forward_values = sorted(v, key=int, reverse=True)
    for k,v in reverse_dict.iteritems():
        reverse_values = sorted(v, key=int, reverse=True)
    results.update({reduced:int(forward_values[-1])+int(reverse_values[-1])})

sorted_dict = sorted(results.iteritems(),key=operator.itemgetter(1))

for item in sorted_dict:
    print item[0]+"\t"+str(item[1])
