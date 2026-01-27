#!/bin/bash 

data=$(dirname $PWD)

[ ! -d $data/Analysis ] && mkdir $data/Analysis
[ ! -d $data/Analysis/0_input_data ] && mkdir -p $data/Analysis/0_input_data

cp $data/4_Picked_OTUs/otu_data/Reotu_table.tsv $data/Analysis/0_input_data
cp $data/4_Picked_OTUs/otu_data/rep-seqs.fasta $data/Analysis/0_input_data
cp $data/4_Picked_OTUs/*rename.qza $data/Analysis/0_input_data
cp $data/6_Taxonomy_gg/taxonomy.qza $data/Analysis/0_input_data/taxonomy_gg.qza
cp $data/6_Taxonomy_silva/taxonomy.qza $data/Analysis/0_input_data/taxonomy_silva.qza
cp $data/6_Taxonomy_gg/taxa_data/taxonomy.tsv $data/Analysis/0_input_data/taxonomy_gg.tsv
cp $data/6_Taxonomy_silva/taxa_data/taxonomy.tsv $data/Analysis/0_input_data/taxonomy_silva.tsv
cp $data/5_Phylogeny_gg/rooted-tree.qza $data/Analysis/0_input_data
cp $data/bin/sample-metadata.txt $data/Analysis/0_input_data

# Rscript $data/bin/A10.1_abundance_data.r $data/Analysis/0_input_data/ $data/Analysis/10_LEfSe/Total/  $data/bin/Total-color.txt 2
