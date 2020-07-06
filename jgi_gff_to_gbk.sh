#!/bin/bash

set -e

# assume .gz files

GFF=$1
FASTA=$2

PREFIX=${GFF%%.*}

echo File prefix: $PREFIX

tmp_dir=$(mktemp -d -t jgi2gbk-XXXXXXXXXX)

zcat $GFF > $tmp_dir/jgi.gff
zcat $FASTA > $tmp_dir/seq.fas

perl ~/Projects/gfftools/jgi2gff.pl $tmp_dir/jgi.gff
~/tmp/genometools-1.6.1/bin/gt gff3 -sort -tidy -addids no $tmp_dir/jgi.gff.gff > $tmp_dir/gt.gff3
~/Projects/gfftools/fix_gendbe_gff.py --no_fix_cds_names --no_fix_escapes_in_attributes --no_include_fasta --gff $tmp_dir/gt.gff3 --fasta $tmp_dir/seq.fas > $tmp_dir/fixed.gff3

~/Projects/gfftools/gff3_to_gbk.py $tmp_dir/fixed.gff3 $tmp_dir/seq.fas
mv $tmp_dir/fixed.gb $PREFIX.gb
rm -rf $tmp_dir
