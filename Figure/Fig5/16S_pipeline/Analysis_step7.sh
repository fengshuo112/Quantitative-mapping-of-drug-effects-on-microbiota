#Plot Alpha Curve
#Current working taxanomic level is Genus

Data=$(dirname $PWD)
Input=$(dirname $PWD)/Analysis/0_input_data/

[ ! -d $Data/Analysis/7_OTUs_alpha ] && mkdir -p $Data/Analysis/7_OTUs_alpha;

Rscript $Data/bin/A7.1_Plot_OTUs_alpha_bar.r $Input $Data/bin/sample-metadata.txt $Data/Analysis/7_OTUs_alpha 4 5
