import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO
import numpy as np

### Sample Command: python LAST_recruitment_plot.py ERR598994_1.fq sample.fna testlastout

inputfile = sys.argv[1] # the FASTQ file to map
reference = sys.argv[2] # the contig file to map against
output = sys.argv[3]    # project prefix for all output files

# CUTOFFS: Change these if you want to make the search more or less stringent
bit_cutoff = 50
alnlength_cutoff = 50
percid_cutoff = 70

# make LAST db
refdb = reference +".lastdb"
cmd = "lastdb "+ refdb +" "+ reference
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

# run LAST search
cmd = "lastal -P 8 "+ refdb +" "+ inputfile +" -f BlastTab -Q 1"
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output+".lastout", "w"))

# get dictionary of contig lengths
contig2length = {}
for record in SeqIO.parse(reference, "fasta"):
	contig2length[record.id] = len(record.seq)

# get added values that should be added on to alignment coordinates
contig_bounds = open(output+".contig_boundaries.txt", "w")
sorted_dict = sorted(contig2length, key=contig2length.get, reverse=True)
added_values = defaultdict(float)
adj = 0
for i in sorted_dict:
	#print(i, contig2length[i])
	added_values[i] = adj
	adj = adj + contig2length[i]
	contig_bounds.write(i +"\t"+ str(added_values[i]) +"\n")
contig_bounds.close()

# parse through lastout file and generate final coordinates and output file
coord_out = open(output+".finalcoords.txt", "w")
coord_out.write("coord\tpercid\n")
handle = open(output+".lastout", "r")
for i in handle.readlines():
	if i.startswith("#"):
		pass
	else:
		line = i.rstrip()
		tabs = line.split("\t")
		contig = tabs[1]
		percid = tabs[2]
		score = tabs[11]
		aln_length = tabs[3]

		if score > bit_cutoff and percid > percid_cutoff and  aln_length > alnlength_cutoff:
			start = float(tabs[8])
			end = float(tabs[9])
			mean = np.mean([start, end])
			final_coord = mean + added_values[contig]
			coord_out.write(str(final_coord) +"\t"+ str(percid) +"\n")
coord_out.close()

# cleanup
os.remove(output+".lastout")
os.remove("stderr.txt")
os.remove("stdout.txt")
for j in os.listdir("."):
	if j.startswith(reference +".lastdb"):
		os.remove(j)












