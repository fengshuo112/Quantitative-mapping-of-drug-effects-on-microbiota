#!/bin/bash

data=$(dirname $PWD)

[ ! -d $data/4_Picked_OTUs ] && mkdir $data/4_Picked_OTUs 

time qiime deblur denoise-16S \
  --i-demultiplexed-seqs $data/1_import_data/demux.qza \
  --p-trim-length 400 \
  --p-sample-stats \
  --o-representative-sequences $data/4_Picked_OTUs/rep-seqs.qza \
  --o-table $data/4_Picked_OTUs/table.qza \
  --o-stats $data/4_Picked_OTUs/deblur-stats.qza

qiime vsearch cluster-features-de-novo \
  --i-table $data/4_Picked_OTUs/table.qza \
  --i-sequences $data/4_Picked_OTUs/rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table $data/4_Picked_OTUs/table-99.qza \
  --o-clustered-sequences $data/4_Picked_OTUs/rep-seqs-99.qza

qiime feature-table summarize \
  --i-table $data/4_Picked_OTUs/table-99.qza \
  --o-visualization $data/4_Picked_OTUs/table-99.qzv \
  --m-sample-metadata-file $data/bin/sample-metadata.txt

time qiime tools export --input-path $data/4_Picked_OTUs/table-99.qzv --output-path $data/4_Picked_OTUs/Diversity_para

python3.8 $data/bin/Diversity_para_cal.py $data/4_Picked_OTUs/Diversity_para/sample-frequency-detail.csv $data/bin/sample-metadata.txt $data/bin/diversity-para.txt
