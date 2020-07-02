#!/usr/bin/env python3

import argparse
import sys
import operator

parser = argparse.ArgumentParser(
  description="""Sort the gff sequences by their length. Assumes that ##sequence-region lines are present""")
parser.add_argument("gff_file", help="The gff3 file")
args = parser.parse_args()

handle = None
if args.gff_file == "-":
  handle = sys.stdin
else:
  handle = open(args.gff_file)

# we store the lines for each sequence as entries in a list: {id: string, length: int, lines: list}
entries = []
seq_id = None
lines = None
sequences_printed = False

def print_sorted_entries():
  for entry in sorted(entries, key=operator.itemgetter("length"), reverse=True):
    print("\n".join(entry['lines']))

for line in handle:
  line = line.strip()
  if line.startswith("##sequence-region"):
    s = line.split()
    seq_len = int(s[3])
    seq_id = s[1]
    entry = {'id': seq_id, 'length': seq_len, 'lines' : [line]}
    lines = entry['lines']
    entries.append(entry)
  elif line.startswith("##FASTA"):
    sequences_printed = True
    print_sorted_entries()
    print(line)
  elif sequences_printed:
    #just print remaining lines
    print(line)
  else:
    if lines != None:
      lines.append(line)

  if not seq_id:
    print(line)

if not sequences_printed:
  print_sorted_entries()
