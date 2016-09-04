#!/bin/bash

if [ $# -lt 2 ]; then
  echo "About:   Create variant sets from annotations"
  echo "Usage:   annot2geneset.sh <SnpSift.jar> <input VCF> <output prefix>"
  exit 1
fi

jar="$1"
vcf="$2"
geneset="$3"

java -Xmx4g -jar $jar extractFields $vcf CHROM POS REF ALT ANN[0].GENE |
  awk -F"\t" 'NR==FNR {x[$1]++} NR>FNR && ($5 in x) {print $1":"$2":"$3":"$4"\t"$5}' $geneset -
