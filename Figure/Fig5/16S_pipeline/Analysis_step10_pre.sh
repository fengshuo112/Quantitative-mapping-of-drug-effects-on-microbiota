#!/bin/bash 

Data=$(dirname $PWD)

[ ! -d $Data/Analysis/10_LEfSe ] && mkdir $Data/Analysis/10_LEfSe;

cp $Data/Analysis/9_PICRUST2/Total.catagrized.kegg.csv $Data/Analysis/0_input_data/Total.kegg.path.csv
Rscript $Data/bin/A10.1_abundance_data.r $Data/Analysis/0_input_data/ $Data/Analysis/10_LEfSe/Total/

for i in $(ls $Data/Analysis/10_LEfSe/*/*.tab);
do
        perl $Data/bin/A9.0_Convert_OTU_table.pl $i;
done

#Prepare the abundance data for LEfSe
for i in $(ls $Data/Analysis/10_LEfSe/*/*.L3.abundance.txt);
do
        ID=$(basename $i);
        dd=$(dirname $i)
        ff=$(basename $i)
        c1=`echo "$ff" | cut -d . -f 1`
        perl $Data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $Data/bin/$c1-color.txt 2;
done

for i in $(ls $Data/Analysis/10_LEfSe/*/*.out.abundance.txt);
do
        ID=$(basename $i);
        dd=$(dirname $i)
        ff=$(basename $i)
        c1=`echo "$ff" | cut -d . -f 1`

        perl $Data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $Data/bin/$c1-color.txt 2;
done

for i in $(ls $Data/Analysis/10_LEfSe/*/*.Genus.abundance.txt);
do
        ID=$(basename $i);
        dd=$(dirname $i)
        ff=$(basename $i)
        c1=`echo "$ff" | cut -d . -f 1`
        perl $Data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $Data/bin/$c1-color.txt 2;
done
