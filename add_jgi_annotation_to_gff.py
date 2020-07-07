#!/usr/bin/env python3
import argparse
import gzip
import sys

parser = argparse.ArgumentParser(description="Adds jgi kegg annotation to feature in a gff file")
parser.add_argument('--gff', help="The gff3 file to correct")
parser.add_argument('--annotation', help="The annotation file (tsv)")
args = parser.parse_args()

annotation_data = {}
with gzip.open(args.annotation, 'rt') as file:
  for line in file:
    if not line.startswith("#"):
      s = line.strip().split('\t')
      id = s[0]
      if id not in annotation_data:
        annotation_data[id] = []
      entry = {'ec': s[1], 'definition': s[2], 'pathway': s[6], 'pathway_class': s[7], 'pathway_type': s[8]}
      annotation_data[id].append(entry)

gff_file = args.gff
# marker that we are in the fasta section (the last section in the file)
in_fasta_section = False


def encodeString(s):
  s = s.replace(",", "%2C").replace("\t", "%09").replace("\n", "%0A").replace("%", "%25").replace(";", "%3B").replace("=", "%3D").replace("&", "%26")
  return s


with open(gff_file) as file:
  for line in file:
    line = line.strip()
    if in_fasta_section:
      print(line)
    else:
      if line.startswith('##FASTA'):
        in_fasta_section = True
        print(line)
      if line.startswith('#'):
        print(line)
      else:
        s = line.split("\t")
        # add product
        if s[2] == "CDS":
          # split attributes
          attrs = s[8].split(';')
          attributes = {}
          for attr in attrs:
            (key, value) = attr.split('=')
            attributes[key] = value

          if 'proteinId' not in attributes:
            print('no proteinId found in:' + line, file=sys.stderr)
          id = attributes['proteinId']
          if id in annotation_data:
            annotation = annotation_data[id]
            text = []
            for entry in annotation:
              text.append(", ".join([entry['ec'], entry['definition'], entry['pathway'], entry['pathway_class'], entry['pathway_class'], entry['pathway_type']]))
            t = "; ".join(text)
            s[8] = s[8] + "; product=" + encodeString(t)
        print("\t".join(s))
