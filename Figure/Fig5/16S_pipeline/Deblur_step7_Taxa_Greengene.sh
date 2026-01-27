#!/bin/bash

data=$(dirname $PWD)

[ ! -d $data/6_Taxonomy_gg ] && mkdir $data/6_Taxonomy_gg
[ ! -d $data/6_Taxonomy_gg/taxa_data ] && mkdir $data/6_Taxonomy_gg/taxa_data
time qiime feature-classifier classify-sklearn \
--i-classifier /home/Metagenome/WZT/16s/greengenes/classifier/gg_13_8_99_classifier.qza \
--i-reads $data/4_Picked_OTUs/rep-seqs-rename.qza \
--o-classification $data/6_Taxonomy_gg/taxonomy.qza

qiime tools export \
--input-path $data/6_Taxonomy_gg/taxonomy.qza \
--output-path $data/6_Taxonomy_gg/taxa_data
