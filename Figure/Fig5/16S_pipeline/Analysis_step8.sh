#!/bin/bash

Data=$(dirname $PWD)
input=$(dirname $PWD)

[ ! -d $Data/Analysis/8_OTUs_beta ] && mkdir -p $Data/Analysis/8_OTUs_beta

qiime tools export \
  --input-path $input/7_Diversity/bray_curtis_distance_matrix.qza \
  --output-path $Data/Analysis/8_OTUs_beta

mv $Data/Analysis/8_OTUs_beta/distance-matrix.tsv $Data/Analysis/8_OTUs_beta/bray-curtis_distance_matrix.tsv

qiime tools export \
  --input-path $input/7_Diversity/jaccard_distance_matrix.qza \
  --output-path $Data/Analysis/8_OTUs_beta
mv $Data/Analysis/8_OTUs_beta/distance-matrix.tsv $Data/Analysis/8_OTUs_beta/jaccard_distance_matrix.tsv

for i in `ls $Data/Analysis/8_OTUs_beta/*distance_matrix.tsv`
do 
	Rscript A7.2_Plot_OTUs_beta_bar.r $i $Data/bin/sample-metadata.txt $Data/Analysis/8_OTUs_beta 4 5
done

