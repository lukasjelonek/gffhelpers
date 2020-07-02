#!/bin/bash
# Calculates the cummulative lengths of features in a gff3 file
# Assumes that each sequence has a ##sequence-region entry with the correct length

GFF_FILE=$1
VERSION=1.0.0
FEATURE_TYPES=("gene" "mRNA" "CDS" "exon" "intron" "intergenic_region")

sep() {
printf "%72s\n" | tr " " "="
}

total_length=$(grep "##sequence-region" $GFF_FILE | awk '{sum += $4}END{print sum}')
echo "Compute gff file statistics v$VERSION"
sep
echo
echo

FORMAT="%39s : %d\n"
printf "%39s : %s\n" "Input file" $GFF_FILE
printf "$FORMAT" "Sequences" $(grep "##sequence-region" $GFF_FILE | wc -l)
printf "%39s : %5sbp (%d bp)\n" "Total length" $(numfmt --to si --format "%2.1f" $total_length) $total_length 
sep

compute_feature_type_total_length () {
    local TYPE=$1
    local length=$(awk -v type=$TYPE '$3==type{sum +=($5-$4+1)}END{print sum}' $GFF_FILE)
    local percentage=$(awk "BEGIN{printf \"%.2f\", $length/$total_length}")
    local CUM_FORMAT="%39s : %5s bp (%9d bp) (%1.2f)\n"
    printf "$CUM_FORMAT" "$TYPE" $(numfmt --to si --format "%2.1f" $length) $length $percentage
}

compute_mean_feature_length () {
    local TYPE=$1
    local length=$(awk -v type=$TYPE '$3==type{sum +=($5-$4+1);count++}END{print sum/count}' $GFF_FILE)
    local MEAN_FORMAT="%39s : %8.2fbp\n"
    printf "$MEAN_FORMAT" "$TYPE" $length
}

compute_feature_counts () {
    local TYPE=$1
    local count=$(awk -v type=$TYPE '$3==type{count++}END{print count}' $GFF_FILE)
    local COUNT_FORMAT="%39s : %8.0f\n"
    printf "$COUNT_FORMAT" "$TYPE" $count
}

compute_complete_cds_mean_length() {
    local length=$(awk '$3 == "CDS"{features[$9]+=$5-$4+1}END{for (x in features ){sum+=features[x]} print sum/length(features)}' $GFF_FILE)
    printf "%39s : %8.2fbp\n" "CDS" $length
}

process_function () {
    local text=$1
    local fnc=$2
    printf "%s\n" "$text"
    for type in ${FEATURE_TYPES[*]}
    do
    $fnc "$type"
    done
    sep
}

process_function "Feature counts" compute_feature_counts
process_function "Cummulative length statistics" compute_feature_type_total_length
process_function "Mean length statistics" compute_mean_feature_length

printf "%s\n" "CDS concatenated mean length"
compute_complete_cds_mean_length
sep
