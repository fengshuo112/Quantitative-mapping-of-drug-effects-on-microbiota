#!/bin/bash

data=$(dirname $PWD)
[ ! -d $data/5_Phylogeny_gg ] && mkdir $data/5_Phylogeny_gg
[ ! -d $data/5_Phylogeny_gg/phyloseq ] && mkdir $data/5_Phylogeny_gg/phyloseq

time qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $data/4_Picked_OTUs/rep-seqs-rename.qza \
  --o-alignment $data/5_Phylogeny_gg/aligned-rep-seqs.qza \
  --o-masked-alignment $data/5_Phylogeny_gg/masked-aligned-rep-seqs.qza \
  --o-tree $data/5_Phylogeny_gg/unrooted-tree.qza \
  --o-rooted-tree $data/5_Phylogeny_gg/rooted-tree.qza

## 导出有根树
qiime tools export \
	--input-path $data/5_Phylogeny_gg/unrooted-tree.qza \
	--output-path $data/5_Phylogeny_gg/phyloseq
mv $data/5_Phylogeny_gg/phyloseq/tree.nwk $data/5_Phylogeny_gg/phyloseq/unrooted_tree.nwk

## 导出无根树
qiime tools export \
--input-path $data/5_Phylogeny_gg/rooted-tree.qza \
--output-path $data/5_Phylogeny_gg/phyloseq
mv $data/5_Phylogeny_gg/phyloseq/tree.nwk $data/5_Phylogeny_gg/phyloseq/rooted_tree.nwk
