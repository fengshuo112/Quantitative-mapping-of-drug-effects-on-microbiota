#!/bin/bash

data=$(dirname $PWD)
out=$(dirname $PWD)
#depth=`head -n 1 $data/bin/diversity-para.txt|cut -d : -f 2`

echo " Calulate  Diversity!"
time qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $data/5_Phylogeny_gg/rooted-tree.qza \
  --i-table $data/4_Picked_OTUs/otu_table-rename.qza \
  --p-sampling-depth 1500 \
  --m-metadata-file $out/bin/sample-metadata.txt \
  --output-dir $out/7_Diversity

echo "Export txt from qza!"
for i in `ls $out/7_Diversity/*.qza`;
do 
	qiime tools export \
		--input-path $i \
		--output-path $out/7_Diversity/diversity_data
done
