#!/usr/bin/env python

"""Find primer conflicts in
a multi-FASTA"""

from __future__ import print_function
from sys import argv
import sys
import subprocess
import optparse
import os

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def parse_blast_file(blast_in,aln_length):
    blast_dict = {}
    for line in open(blast_in, "rU"):
        newline = line.strip()
        fields = newline.split()
        if int(fields[3])>int(aln_length) and int(fields[5])==0:
            try:
                blast_dict[fields[0]].append(fields[1])
            except KeyError:
                blast_dict[fields[0]] = [fields[1]]
    print("Problematic primers with an alignment length of %s, number of bad alignments" % aln_length)
    for k,v in blast_dict.iteritems():
        if len(v)>1:
            #print(k,len(v)-1)
            print(k,"\t".join(v))
def main(primer_fasta,aln_length):
    dependencies = ['blastn','makeblastdb']
    for dependency in dependencies:
        ra = subprocess.call(['which', '%s' % dependency])
        if ra == 0:
            pass
        else:
            print("%s is not in your path, but needs to be!" % dependency)
            sys.exit()
    primer_path = os.path.abspath(primer_fasta)
    subprocess.check_call("makeblastdb -in %s -dbtype nucl > /dev/null 2>&1" % primer_path, shell=True)
    subprocess.check_call("blastn -task blastn -query %s -db %s -perc_identity 100 -dust no -num_threads 4 -evalue 20 -word_size 4 -outfmt 6 -out self_blast.out > /dev/null 2>&1" % (primer_path,primer_path), shell=True)
    parse_blast_file("self_blast.out",aln_length)
    #os.system("rm self_blast.out")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-f", "--primer_fasta", dest="primer_fasta",
                      help="path to primers in multi-FASTA [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-l", "--aln_length", dest="aln_length",
                      help="minimum alignment length, defaults to 10",
                      type="int", action="store", default="10")
    options, args = parser.parse_args()
    mandatories = ["primer_fasta"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.primer_fasta,options.aln_length)
