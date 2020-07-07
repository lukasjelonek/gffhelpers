# gffhelpers
Helper scripts to handle genome annotation in gff format

This Repository contains scripts that helped me handling gff files. All these
scripts had specific purpose and cannot be generalized. So use them with
caution and adapt them to your needs.

## Tools

* `add_intergenic_regions.py` - Adds intergenic regions between 'sequence
  start':'first gene', non-overlapping genes and 'last gene':'sequence end'
* `add_jgi_annotation_to_gff.py` - Includes jgi Kegg.tsv annotation with a gff3
  file, based on the proteinId attribute
* `compute_feature_density.sh` - Computes statistics in a gff3 file: sequence
  count, feature count per type, mean/median feature length, cummulative
  feature length, aggregated cds mean/median length
* `compute_intron_histogram_from_gff.sh` - 
* `fix_gendbe_gff.py` - Processes gendbe-gff3 output. Escapes unescaped
  characters in attributes, removes numbers from cds ids, adds fasta sequence
  to gff, adds sequence regions based on fasta
* `gff3_to_gbk.py` - Converts a gff3 file to a genbank file.
* `jgi2gff3.py` - Converts a jgi gff file (version 1, 2 or gtf?) to gff3
* `jgi_gff_to_gbk.sh` - A pipeline to convert jgi data to genbank format
* `sort_gff_by_sequence_region_length.py` - Sorts the entries in a gff file by
  the sequence length (longest first)

## Dependencies

* `bash`
* `awk`
* `sed`
* `python`
* `biopython`
* `genometools`
