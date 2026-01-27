#!/bin/bash

data=$(dirname $PWD)
Taxa=/home/Metagenome/Ref_database/SILVA/classifiers

[ ! -d $data/6_Taxonomy_silva ] && mkdir $data/6_Taxonomy_silva
[ ! -d $data/6_Taxonomy_silva/taxa_data ] && mkdir $data/6_Taxonomy_silva/taxa_data
time qiime feature-classifier classify-sklearn \
	--i-classifier $Taxa/silva_138_99_nb_classifier.qza \
	--i-reads $data/4_Picked_OTUs/rep-seqs-rename.qza \
	--o-classification $data/6_Taxonomy_silva/taxonomy.qza

qiime tools export \
	--input-path $data/6_Taxonomy_silva/taxonomy.qza \
	--output-path $data/6_Taxonomy_silva/taxa_data
