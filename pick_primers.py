#!/usr/bin/env python

"""Picks primers from a given gene,
making sure that it hits all target genomes
and misses as many non-target genomes
as possible"""

from __future__ import print_function
import optparse
import sys
import os
import subprocess
import glob
import collections
from collections import OrderedDict
try:
    from Bio import SeqIO
except:
    print("BioPython is not installed, but needs to be")
    sys.exit()

"""Tests for file inputs"""
def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print('%s file cannot be opened' % option)
        sys.exit()

def test_dir(option, opt_str, value, parser):
    if os.path.exists(value):
        setattr(parser.values, option.dest, value)
    else:
        print("directory of fastas cannot be found")
        sys.exit()
"""Finished with tests for file inputs"""

def get_seq_name(fasta_in):
    name = os.path.basename(fasta_in)
    return name

def combine_seqs(dir_path):
    handle = open("combined.seqs", "w")
    for infile in glob.glob(os.path.join(dir_path, '*.fasta')):
        names = get_seq_name(infile)
        reduced = names.replace('.fasta','')
        handle.write(">"+str(reduced)+"\n")
        for record in SeqIO.parse(open(infile), "fasta"):
            print(record.seq, file=handle)
    handle.close()

def modify_input_file(config_file, gene, gene_name, gene_length, upper, lower):
    infile = open(gene, "U")
    outfile = open("sequence.txt", "w")
    sequence = []
    gene_min = int(lower)
    gene_max = int(upper)
    for record in SeqIO.parse(infile, "fasta"):
        if record.id == gene_name:
            outfile.write(str(record.seq))
            marker_name = record.id
    infile.close()
    outfile.close()
    for line in open("sequence.txt", "U"):
        sequence.append(line)
    config_in = open(config_file, "U")
    config_out = open("config_modified.txt", "w")
    for line in config_in:
        if line.startswith("SEQUENCE_ID"):
            config_out.write("SEQUENCE_ID=%s" % marker_name+"\n")
        elif line.startswith("SEQUENCE_TEMPLATE"):
            config_out.write("SEQUENCE_TEMPLATE=%s" % "".join(sequence)+"\n")
        elif line.startswith("PRIMER_NUM_RETURN"):
            fields=line.split("=")
            num_to_return = fields[1]
	elif line.startswith("PRIMER_PRODUCT_SIZE_RANGE"):
	    config_out.write("PRIMER_PRODUCT_SIZE_RANGE=%s-%s" % (gene_min,gene_max)+"\n")
        else:
            config_out.write(line)
    config_in.close()
    config_out.close()
    return num_to_return

def process_primer3_output(primer3_file, num_records):
    outfile = open("all_primers.fasta", "w")
    for i in range(int(num_records)):
        for line in open(primer3_file, "U"):
            if line.startswith("PRIMER_LEFT_%s=" % i):
                temp_fields = line.split("=")
                name_fields = temp_fields[1].split(",")
                forward_name = name_fields[0]
            elif line.startswith("PRIMER_RIGHT_%s=" % i):
                temp_fields = line.split("=")
                name_fields = temp_fields[1].split(",")
                reverse_name = name_fields[0]
        for line in open(primer3_file, "U"):
             if line.startswith("PRIMER_LEFT_%s_SEQUENCE" % i):
                fields = line.split("=")
                outfile.write(">%s_forward" % forward_name+"\n")
                outfile.write(fields[1])
             elif line.startswith("PRIMER_RIGHT_%s_SEQUENCE" % i):
                fields = line.split("=")
                outfile.write(">%s_reverse" % reverse_name+"\n")
                outfile.write(fields[1])
    outfile.close()

def create_separate_databases(combined_in, target_ids):
    infile = open(combined_in, "U")
    data = open(target_ids, "U").read().splitlines()
    num_targets = len(data)
    output_handle = open("target_database.seqs", "w")
    output2_handle = open("non_target_database.seqs", "w")
    targets=[ ]
    non_targets=[ ]
    for record in SeqIO.parse(infile, "fasta"):
        if record.id in data:
            targets.append(record)
        else:
            non_targets.append(record)
    SeqIO.write(targets, output_handle, "fasta")
    SeqIO.write(non_targets, output2_handle, "fasta")
    infile.close()
    output_handle.close()
    output2_handle.close()
    return num_targets

def parse_blast_report(infile, num_targets, primer_dict):
    my_dict = {}
    outfile = open("target_hit_results.txt", "w")
    for line in open(infile, "U"):
        fields = line.split()
        """Checks that all primers match 100 percent and also are the same length as the alignment"""
        if fields:
            if int(primer_dict.get(fields[0])) == int(fields[3]) and float(fields[2]) == 100:
                try:
                    my_dict[fields[0]].append(fields[3])
                except KeyError:
                    my_dict[fields[0]]=[fields[3]]
    	else:
    	    print("no hits found")
    for k,v in my_dict.iteritems():
        my_values = sorted(v, key=int, reverse=True)
        try:
            counter=collections.Counter(my_values[0:num_targets])
        except:
	        print(k,my_values[0],"0",file=outfile)
        values=counter.values()
        """this makes sure that all hits have the same alignment length
        should now be redundant, but shouldn't hurt anything?"""
        if len(my_values)==num_targets:
            print(k,my_values[0],my_values[num_targets-1],file=outfile)
        elif len(my_values)>num_targets:
	        print(k,my_values[0],my_values[num_targets],file=outfile)
        else:
	        pass
    outfile.close()

def parse_non_target_blast(infile, primer_names):
    my_dict = {}
    outfile = open("non_target_hit_results.txt", "w")
    for line in open(infile, "U"):
        fields = line.split()
        try:
            my_dict[fields[0]].append(fields[3])
        except KeyError:
            my_dict[fields[0]]=[fields[3]]
    for primer_name in primer_names:
        if primer_name in my_dict:
            pass
        else:
            my_dict.update({primer_name:"0"})
    for k,v in my_dict.iteritems():
        my_values = sorted(v, key=int, reverse=True)
        print(k,my_values[0],file=outfile)
    outfile.close()

def derep_primers(infile, outfile):
    records = [ ]
    record_ids = []
    outfile = open(outfile, "w")
    for record in SeqIO.parse(open(infile, "U"), "fasta"):
        if record.id in record_ids:
            pass
        else:
            records.append(record)
            record_ids.append(record.id)
    SeqIO.write(records, outfile, "fasta")
    outfile.close()
    return record_ids

def report_results(target, non_target, primer_names, gene_name):
    target_dict = {}
    non_target_dict = {}
    outfile = open("%s_results.txt" % gene_name, "w")
    for line in open(target, "U"):
        fields = line.split()
        target_dict.update({fields[0]:fields[2]})
    for line in open(non_target, "U"):
        fields = line.split()
        non_target_dict.update({fields[0]:fields[1]})
    sorted_targets = list(sorted(target_dict, key=target_dict.__getitem__))
    sorted_non_targets = list(sorted(non_target_dict, key=non_target_dict.__getitem__))
    target_results_dict={}
    non_target_results_dict={}
    results = {}
    non_target_totals={}
    for primer_name in primer_names:
        if primer_name in target_dict:
	    target_results_dict.update({primer_name:sorted_targets.index(primer_name)})
	else:
	    pass
	if primer_name in non_target_dict and primer_name in target_dict:
            non_target_results_dict.update({primer_name:sorted_non_targets.index(primer_name)})
            non_target_totals.update({primer_name:int(target_dict.get(primer_name))+int(non_target_dict.get(primer_name))})
	elif primer_name in target_dict and primer_name not in non_target_dict:
	    non_target_results_dict.update({primer_name:"0"})
            non_target_totals.update({primer_name:int(target_dict.get(primer_name))+int("0")})
	else:
	    pass
    for primer_name in primer_names:
        if primer_name in target_dict:
            results.update({primer_name:int(target_results_dict.get(primer_name))+int(non_target_results_dict.get(primer_name))})
	else:
	    pass
    sorted_results = [(key, str(val)) for val, key in sorted((int(val), key) for key, val in results.iteritems())]
    for result in sorted_results:
        print("\t".join(result)+"\t"+str(non_target_totals.get(result[0])),file=outfile)

def get_gene_name(gene_path):
    for record in SeqIO.parse(open(gene_path, "U"), "fasta"):
        return record.id, len(record.seq)

def get_targets(ids):
    data = open(ids, "U").read().splitlines()
    num_targets = len(data)
    return num_targets

def get_primer_lengths(primers_in):
    primer_dict = {}
    for record in SeqIO.parse(open(primers_in, "rU"), "fasta"):
        primer_dict.update({record.id:len(record.seq)})
    return primer_dict

def main(config_file, gene, target_ids, directory, upper, lower):
    """test for dependencies"""
    dependencies = ['primer3_core','blastall','makeblastdb']
    for dependency in dependencies:
        ra = subprocess.call(['which', '%s' % dependency])
        if ra == 0:
            pass
        else:
            print("%s is not in your path, but needs to be!" % dependency)
            sys.exit()
    """dependency testing finished"""
    start_dir = os.getcwd()
    start_path = os.path.abspath("%s" % start_dir)
    config_file_path = os.path.abspath("%s" % config_file)
    dir_path = os.path.abspath("%s" % directory)
    gene_path = os.path.abspath("%s" % gene)
    target_ids = os.path.abspath("%s" % target_ids)
    """start pipeline now"""
    for record in SeqIO.parse(open(gene_path), "fasta"):
        gene_name = record.id
        gene_length = len(record.seq)
        if int(gene_length)<int(upper):
            pass
        else:
            num_records = modify_input_file(config_file, gene_path, gene_name, gene_length, upper, lower)
            subprocess.check_call("primer3_core -output=primer3_out < config_modified.txt", shell=True)
            process_primer3_output("primer3_out", num_records)
            if os.path.isfile("combined.seqs"):
                pass
            else:
                combine_seqs(directory)
            if os.path.isfile("target_database.seqs"):
                num_targets = get_targets(target_ids)
                pass
            else:
                num_targets = create_separate_databases("combined.seqs", target_ids)
            os.system("makeblastdb -in target_database.seqs -dbtype nucl > /dev/null 2>&1")
            os.system("makeblastdb -in non_target_database.seqs -dbtype nucl > /dev/null 2>&1")
            primer_names = derep_primers("all_primers.fasta", "reduced_primers.fasta")
            primer_dict = get_primer_lengths("reduced_primers.fasta")
            #os.system("blastall -p blastn -i reduced_primers.fasta -d target_database.seqs -b 2000 -v 2000 -e 10 -m 8 -o target_blast.out")
            os.system("blastn -task blastn -query reduced_primers.fasta -db target_database.seqs -num_alignments 2000 -outfmt 6 -out target_blast.out -evalue 10")
            """check this function"""
            parse_blast_report("target_blast.out", num_targets, primer_dict)
            #os.system("blastall -p blastn -i reduced_primers.fasta -d non_target_database.seqs -b 2000 -v 2000 -e 10 -m 8 -o non_target_blast.out")
            os.system("blastn -task blastn -query reduced_primers.fasta -db non_target_database.seqs -num_alignments 2000 -outfmt 6 -out non_target_blast.out -evalue 10")
            parse_non_target_blast("non_target_blast.out", primer_names)
            report_results("target_hit_results.txt", "non_target_hit_results.txt", primer_names, gene_name)
            os.system("cp reduced_primers.fasta %s_primers.seqs" % gene_name)
            os.system("rm reduced_primers.fasta target_blast.out target_hit_results.txt non_target_blast.out non_target_hit_results.txt all_primers.fasta *.for *.rev config_modified.txt")

if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-c", "--config_file", dest="config_file",
                      help="/path/to/primer3 config file [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-g", "--gene", dest="gene",
                      help="/path/to/gene sequence file [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-t", "--target", dest="target_ids",
                      help="/path/to/text file with target ids [REQUIRED]",
                      type="string", action="callback", callback=test_file)
    parser.add_option("-d", "--directory", dest="directory",
                      help="/path/to/directory of genomes in FASTA format [REQUIRED]",
                      type="string", action="callback", callback=test_dir)
    parser.add_option("-u", "--upper", dest="upper",
                      help="upper range for amplicon size, defaults to 400",
                      type="int", action="store", default="400")
    parser.add_option("-l", "--lower", dest="lower",
                      help="lower range for amplicon size, defaults to 200",
                      type="int", action="store", default="200")

    options, args = parser.parse_args()
    mandatories = ["config_file", "gene", "target_ids", "directory"]
    for m in mandatories:
        if not getattr(options, m, None):
            print("\nMust provide %s.\n" %m)
            parser.print_help()
            exit(-1)

    main(options.config_file,options.gene,options.target_ids,options.directory,options.upper,options.lower)
