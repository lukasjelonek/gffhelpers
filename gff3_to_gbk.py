#!/usr/bin/env python3
"""Convert a GFF and associated FASTA file into GenBank format. Only includes genes, mrnas and joined cds

Usage:
gff3_to_gbk.py <GFF annotation file> <FASTA sequence file>
"""
import sys
import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import CompoundLocation, SeqFeature

def main(gff_file, fasta_file):
    out_file = "%s.gb" % os.path.splitext(gff_file)[0]
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    records = []
    for rec in gff_iter:
      features = []
      for f in rec.features:
        features.append(f)
        if f.type == "gene":
          for mrna in _filter_subfeatures(f, 'mRNA'):
            exons = _join_features(_filter_subfeatures(mrna, 'CDS'))
            mrna.location = exons.location
            features.append(mrna)
            features.append(_join_features(_filter_subfeatures(mrna, 'CDS')))
            features.extend(_filter_subfeatures(mrna, 'exon'))
      rec.features = features
      records.append(rec)
    SeqIO.write(records, out_file, "genbank")

def _filter_subfeatures(feature, type):
  return list(filter(lambda f: f.type == type, feature.sub_features))

def _join_features(features): 
  if len(features) == 1:
    return features[0]
  else:
    locations = list(map(lambda f: f.location, features))
    loc = CompoundLocation(locations)
    return SeqFeature(
      location=loc, 
      type=features[0].type, 
      strand=features[0].strand,
      id=features[0].id,
      qualifiers=features[0].qualifiers,
      ref=features[0].ref,
      ref_db=features[0].ref_db)

if __name__ == "__main__":
    main(*sys.argv[1:])
