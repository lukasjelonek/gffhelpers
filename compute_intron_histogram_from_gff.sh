#!/bin/bash
# Generates a histogram of introns per gene
#
# assumptions
#  * gff contains intron features
#  * each intron only has a Parent entry in its additional attributes
#  * all introns to one gene are next to each other


GFF_FILE=$1
awk '$3=="intron" {print $9}' $GFF_FILE \
    | uniq -c \
    | awk '{print $1}' \
    | sort -n \
    | less \
    | uniq -c \
    | awk 'BEGIN{OFS="\t";print "Introns", "count"}{print $2,$1}'
