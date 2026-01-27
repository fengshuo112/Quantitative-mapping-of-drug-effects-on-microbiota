#!/bin/bash

data=$(dirname $PWD)

[ ! -d $data/4_Picked_OTUs/otu_data ] && mkdir $data/4_Picked_OTUs/otu_data
# 导出特征表
qiime tools export --input-path $data/4_Picked_OTUs/uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza --output-path $data/4_Picked_OTUs/otu_data
less $data/4_Picked_OTUs/otu_data/dna-sequences.fasta |paste - -|sed '1i ASVID,seq' > $data/4_Picked_OTUs/otu_data/rep.fa

Rscript ModifyRepSeqname.r $data/4_Picked_OTUs/otu_data

less $data/4_Picked_OTUs/otu_data/rep.xls|sed 's/\r//g'|tr "\t" "\n" > $data/4_Picked_OTUs/otu_data/rep-seqs.fasta

# 导入Qiime2中
time qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $data/4_Picked_OTUs/otu_data/rep-seqs.fasta \
--output-path $data/4_Picked_OTUs/rep-seqs-rename.qza
