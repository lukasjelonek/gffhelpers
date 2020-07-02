#!/usr/bin/env python3
from Bio import SeqIO
import re
import argparse

parser = argparse.ArgumentParser(description="Helper script to fix some errors in gendb-e generated gff files")
parser.add_argument('--gff', help="The gff3 file to correct")
parser.add_argument('--fasta', help="The fasta file. Is used to retrieve the contig lengths")
parser.add_argument('--no_add_sequence_regions', help="Add '##sequence_region'-lines for each sequence. Requires the fasta file.")
parser.add_argument('--no_include_fasta', help="Include the fasta sequence in the gff.")
parser.add_argument('--no_fix_escapes_in_attributes', help="If possible escapes '=,&' in the attribute-values.")
parser.add_argument('--no_fix_cds_names', help="removes trailing numbers in the cds identifiers. Only useful for gendb-e output.")
args = parser.parse_args()

fasta_file = args.fasta

add_sequence_region_entries = not args.no_add_sequence_regions
fix_escape_attributes = not args.no_fix_escapes_in_attributes
fix_cds_names = not args.no_fix_cds_names
include_fasta = not args.no_include_fasta

contig_data = {}
if add_sequence_region_entries:
  with open(fasta_file) as file:
    for record in SeqIO.parse(file, "fasta"):
      contig_data[record.id] = len(record)

gff_file = args.gff
print("##gff-version 3")

with open(gff_file) as file:
  last_id = ""
  for line in file:
    if not line.startswith("#"):
      s = line.strip().split("\t", 9)
      cur_id = s[0]

      ## add sequence-regions
      if add_sequence_region_entries:
        if (cur_id != last_id):
          if cur_id in contig_data:
            print("##sequence-region " + cur_id + " 1 " + str(contig_data[cur_id]))
          else:
            raise Exception("Sequence '" + cur_id + "' not found in fasta")

      last_id = cur_id
      
      cur_type = s[2]
      if fix_escape_attributes:
        if cur_type == "mRNA":
          attributes = s[8].split(';')
          corrected_attributes = [] 
          for attribute in attributes:
            (key, value) = attribute.split("=", 1)
            value = value.replace("=", "%3D").replace(",", "%2C").replace("&", "%26")
            corrected_attributes.append(key + "=" +value)
          s[8] = ";".join(corrected_attributes)
        
      if fix_escape_attributes:
        if s[2] == 'CDS':
          pattern = re.compile('ID=(.*\.cds)\d+;?')
          s[8] = pattern.sub(r'ID=\1;',s[8])
        
      print("\t".join(s))


if include_fasta:
  with open(fasta_file) as file:
    print('##FASTA')
    for line in file:
      print(line.strip())
