#!/bin/bash

data=$(dirname $PWD)
input=$(dirname $PWD)
[ ! -d $data/Analysis ] && mkdir  $data/Analysis
[ ! -d $data/Analysis/1_Abundance ] && mkdir -p $data/Analysis/1_Abundance
[ ! -d $data/Analysis/2_PcoA ] && mkdir -p $data/Analysis/2_PcoA
arg1=("B_ctl-M_ctl" "B_ctl-M_D" "B_ctl-M_E" "M_ctl-M_D" "M_ctl-M_E" "M_D-M_E")

Rscript qiime2_barplot_abundance.r $input/Analysis/0_input_data/ "B_ctl-M_ctl-M_D-M_E" $data/Analysis/1_Abundance/
Rscript A1.2.2_Plot_Abundance_between.r $input/Analysis/0_input_data/ $data/Analysis/1_Abundance/ $data/bin/sample-metadata.txt 4 15 $arg1
Rscript A1.3_Plot_Pcoa.r $input/Analysis/0_input_data/ $data/Analysis/2_PcoA/
