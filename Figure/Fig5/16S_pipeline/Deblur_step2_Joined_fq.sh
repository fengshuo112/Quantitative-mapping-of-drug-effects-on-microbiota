#!/bin/bash

data=$(dirname $PWD)

[ ! -d $data/2_combined_fq ] && mkdir $data/2_combined_fq

time qiime vsearch join-pairs \
  --i-demultiplexed-seqs $data/1_import_data/demux.qza \
  --o-joined-sequences $data/2_combined_fq/demux-joined.qza

time qiime tools export --input-path $data/2_combined_fq/demux-joined.qza --output-path $data/2_combined_fq

gunzip $data/2_combined_fq/*.gz
