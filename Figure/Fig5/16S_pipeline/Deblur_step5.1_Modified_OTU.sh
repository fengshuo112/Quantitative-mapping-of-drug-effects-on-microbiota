#!/bin/bash

# 导出特征表并修改名称
# 路线：# otu_table.qza --> feature-table.biom ---> otu_table.tsv

data=$(dirname $PWD)
[ ! -d $data/4_Picked_OTUs ] && mkdir $data/4_Picked_OTUs
[ ! -d $data/4_Picked_OTUs/otu_data ] && mkdir $data/4_Picked_OTUs/otu_data
qiime tools export --input-path $data/4_Picked_OTUs/uchime-dn-out/table-nonchimeric-wo-borderline.qza --output-path $data/4_Picked_OTUs/otu_data

biom convert -i $data/4_Picked_OTUs/otu_data/feature-table.biom -o $data/4_Picked_OTUs/otu_data/otu_table.tsv --to-tsv
sed -i '1d' $data/4_Picked_OTUs/otu_data/otu_table.tsv
sed -i 's/#OTU ID/ASV/' $data/4_Picked_OTUs/otu_data/otu_table.tsv

Rscript ModifyOTUname.r $data/4_Picked_OTUs/otu_data

# 导入特征表
# 路线：# otu_table.tsv --> feature-table.biom ---> otu_table.qza
biom convert -i $data/4_Picked_OTUs/otu_data/Reotu_table.tsv -o $data/4_Picked_OTUs/otu_data/feature-table.biom --to-hdf5 --table-type="OTU table"

time qiime tools import \
  --input-path $data/4_Picked_OTUs/otu_data/feature-table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path $data/4_Picked_OTUs/otu_table-rename.qza
