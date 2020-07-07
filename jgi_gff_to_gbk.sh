#!/bin/bash

set -e

# assume .gz files

GFF=$1
FASTA=$2
KEGG_ANNOTATION=$3
ORGANISM_NAME=$4
STRAIN=$5

PREFIX=${GFF##*/}
PREFIX=${PREFIX%%.*}

echo File prefix: $PREFIX

tmp_dir=$(mktemp -d -t jgi2gbk-XXXXXXXXXX)

function finish {
    echo Cleanup
    rm -rf $tmp_dir
}

trap finish EXIT

echo Extracting fasta
zcat $FASTA > $tmp_dir/seq.fas

echo Convert to gff3
~/Projects/gfftools/jgi2gff3.py $GFF > $tmp_dir/jgi.gff
echo Sort and cleanup gff3
~/tmp/genometools-1.6.1/bin/gt gff3 -sort -tidy -addids no -retainids $tmp_dir/jgi.gff > $tmp_dir/gt.gff3
echo Add sequence-regions
~/Projects/gfftools/fix_gendbe_gff.py --no_fix_cds_names --no_fix_escapes_in_attributes --no_include_fasta --gff $tmp_dir/gt.gff3 --fasta $tmp_dir/seq.fas > $tmp_dir/fixed.gff3
echo Add Kegg annotation
~/Projects/gfftools/add_jgi_annotation_to_gff.py --gff $tmp_dir/fixed.gff3 --annotation $KEGG_ANNOTATION > $tmp_dir/annotation.gff3

echo Convert to gbk
~/Projects/gfftools/gff3_to_gbk.py $tmp_dir/annotation.gff3 $tmp_dir/seq.fas "$ORGANISM_NAME" "$STRAIN"

mv $tmp_dir/annotation.gb $PREFIX.gb
mv $tmp_dir/annotation.gff3 $PREFIX.gff3

