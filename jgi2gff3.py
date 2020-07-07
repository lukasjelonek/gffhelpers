#!/usr/bin/env python3
import gzip
import sys

# converts a jgi gff file to gff3
# creates mrna and gene features with the names of the cds
# removes start and stop codons as they are included in the cds and ommitted in gff3
# usage: jgi2gff3.py <gff-file>
# output: gff3 to stdout

def read_entry(line):
  line = line.strip()
  s = line.split("\t")
  entry = {
    'seqid': s[0],
    'source': s[1],
    'type': s[2],
    'start': int(s[3]),
    'end': int(s[4]),
    'score': s[5],
    'strand': s[6],
    'phase': s[7]
  }
  attr = s[8].split("; ")
  attributes = {}
  for a in attr:
    (key, value) = a.split(" ", 1)
    attributes[key] = value
  entry['attributes'] = attributes
  return entry


def _filter_by_types(entries, types):
  'return the entries that have one of the given types'
  return list(filter(lambda e: any(map(lambda t: e['type'] == t, types)), entries))


def _remove_by_types(entries, types):
  'return the entries that do not have one of the given types'
  return list(filter(lambda e: not any(map(lambda t: e['type'] == t, types)), entries))


def _print_entry(entry):
  print("\t".join([
    entry['seqid'],
    entry['source'],
    entry['type'],
    str(entry['start']),
    str(entry['end']),
    entry['score'],
    entry['strand'],
    entry['phase'],
    ";".join(map(lambda x: x[0] + "=" + x[1], entry['attributes'].items()))
  ]))


def _check_that_transcript_ids_are_the_same(entries):
  ids = set(map(lambda x: x['attributes']['transcriptId'], _filter_by_types(entries, ['exon'])))
  if len(ids) != 1:
    raise Exception("transcriptIds differ: " + str(ids))


def _check_that_protein_ids_are_the_same(entries):
  ids = set(map(lambda x: x['attributes']['proteinId'], _filter_by_types(entries, ['CDS'])))
  if len(ids) != 1:
    raise Exception("proteinIds differ: " + str(ids))


def process_entry(entries):
  _check_that_transcript_ids_are_the_same(entries)
  _check_that_protein_ids_are_the_same(entries)
  name = entries[0]['attributes']['name'].replace('"', '')
  # find min start and max end for mRNA and gene feature
  start = min(map(lambda e: e['start'], entries))
  end = max(map(lambda e: e['end'], entries))
  # create mRNA
  mRNA_entry = {
    'seqid': entries[0]['seqid'],
    'source': entries[0]['source'],
    'type': 'mRNA',
    'start': start,
    'end': end,
    'score': entries[0]['score'],
    'strand': entries[0]['strand'],
    'phase': '.',
    'attributes': {'ID': name+".t", 'Parent': name+'.g'}
  }
  # create gene
  gene_entry = {
    'seqid': entries[0]['seqid'],
    'source': entries[0]['source'],
    'type': 'gene',
    'start': start,
    'end': end,
    'score': entries[0]['score'],
    'strand': entries[0]['strand'],
    'phase': '.',
    'attributes': {'ID': name+".g"}
  }

  # fix attributes for all entries
  for e in entries:
    e['attributes']['Parent'] = name+".t"
    e['attributes'].pop('name')

  # drop start and stop codons
  entries = _remove_by_types(entries, ['start_codon', 'stop_codon'])

  _print_entry(gene_entry)
  _print_entry(mRNA_entry)
  for e in entries:
    _print_entry(e)

  print('###')


filename = sys.argv[1]
processed_entries = {}
entries = None
last_name = None
print("##gff-version 3")
with gzip.open(filename, "rt") as file:
  for line in file:
    entry = read_entry(line)
    if entry['attributes']['name'] != last_name:
      if last_name:
        process_entry(entries)
        if last_name in processed_entries:
          raise Exception("Multiple entries for " + last_name)
        processed_entries[last_name] = 1
      last_name = entry['attributes']['name']
      entries = [entry]
    else:
      entries.append(entry)
