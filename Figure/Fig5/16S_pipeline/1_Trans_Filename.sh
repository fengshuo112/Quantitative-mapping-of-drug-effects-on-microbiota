#!/bin/bash

data=$(dirname $PWD)
#gunzip $data/*.gz
cd $data/0_raw_data
: > Trans.sh

ls *R1_001.fastq|awk -F "_" '{print "mv "$0" "$1"_16s_R1.fastq"}' >> Trans.sh
ls *R2_001.fastq|awk -F "_" '{print "mv "$0" "$1"_16s_R2.fastq"}' >> Trans.sh

bash Trans.sh
rm Trans.sh

