#!/usr/bin/env python3

import argparse
# computes intergenic_regions in gff3 files
# assumptions:
#  * gff contains ##sequence-region lines
#  * gff contains gene lines
#  * gene lines are sorted by seq_id and start

parser = argparse.ArgumentParser(
  description="""Compute intergenic_regions in a gff3 file. 
  Assumes that the gff file contains '##sequence-region' lines and that the genes are sorted by sequence and start-coordinate.
  Will create a new 'intergenic_region'-entry between contig-start and first gene, between non-overlapping genes and between last-gene and contig-end""", 
  epilog="The resulting gff3 is written to stdout.")
parser.add_argument("gff_file", help="The gff3 file")
args = parser.parse_args()

start = None
seq_id = None
seq_len = None
count = None

def print_intergenic_region(id, start, end, number):
  region = [id, 'add_intergenic_regions.py', 'intergenic_region', str(start), str(end), ".", ".", ".", "ID="+id+".intergenic." + str(number)]
  print("\t".join(region))

with open(args.gff_file) as file:
  for line in file:
    if line.startswith("##sequence-region"):
      ## add a intergenic_region at the end of the contig
      if seq_id and start < seq_len:
        print_intergenic_region(seq_id, start, seq_len, count)

      # reset the variables for the next contig
      s = line.split()
      seq_len = int(s[3])
      start = 1
      seq_id = s[1]
      count = 1

    if line.startswith("#"):
      # just print the comment lines
      print(line.strip())
    else:
      s = line.split("\t", 9)
      if s[2] == 'gene':
        # add a intergenic_region befor each gene that does not overlap the previous one
        gene_start = int(s[3])
        gene_end = int(s[4])
        if start < gene_start-1:
          print_intergenic_region(seq_id, start, gene_start-1, count)
          count = count + 1
        start = gene_end + 1 
      # print each feature line
      print(line.strip())

  # don't miss the last intergenic region
  if seq_id and start < seq_len:
    print_intergenic_region(seq_id, start, seq_len, count)
