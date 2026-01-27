#!/bin/bash

data=$(dirname $PWD)

[ ! -d $data/3_filter_fq ]  && mkdir $data/3_filter_fq
[ ! -d $data/3_filter_fq/filter_data ] && mkdir -p $data/3_filter_fq/filter_data

time qiime quality-filter q-score \
  --i-demux $data/2_combined_fq/demux-joined.qza \
  --o-filtered-sequences $data/3_filter_fq/demux-joined-filtered.qza \
  --o-filter-stats $data/3_filter_fq/demux-joined-filter-stats.qza

time qiime tools export --input-path $data/3_filter_fq/demux-joined-filtered.qza --output-path $data/3_filter_fq
gunzip $data/3_filter_fq/*.gz
for i in `ls $data/3_filter_fq/*.fastq`;
do
	seqtk seq -A $i > $data/3_filter_fq/$(basename $i .fastq).fa
done
rm -rf $data/3_filter_fq/filter_data
