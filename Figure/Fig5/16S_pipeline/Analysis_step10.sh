#!/bin/bash
data=$(dirname $PWD)

[ ! -d $data/Analysis/10_LEfSe/Total ] && mkdir -p $data/Analysis/10_LEfSe/Total;

cp $data/Analysis/9_PICRUST2/Total.catagrized.kegg.csv $data/Analysis/0_input_data/Total.kegg.path.csv
Rscript $data/bin/A10.1_abundance_data.r $data/Analysis/0_input_data/ $data/Analysis/10_LEfSe/Total/

for i in $(ls $data/Analysis/10_LEfSe/*/*.tab);
do
        perl $data/bin/A9.0_Convert_OTU_table.pl $i;
done

#Prepare the abundance data for LEfSe
for i in $(ls $data/Analysis/10_LEfSe/*/*.L3.abundance.txt);
do
        ID=$(basename $i);
        dd=$(dirname $i)
        ff=$(basename $i)
        c1=`echo "$ff" | cut -d . -f 1`
        perl $data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $data/bin/$c1-color.txt 2;
done

for i in $(ls $data/Analysis/10_LEfSe/*/*.out.abundance.txt);
do
        ID=$(basename $i);
        dd=$(dirname $i)
        ff=$(basename $i)
        c1=`echo "$ff" | cut -d . -f 1`

        perl $data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $data/bin/$c1-color.txt 2;
done

for i in $(ls $data/Analysis/10_LEfSe/*/*.Genus.abundance.txt);
do
        ID=$(basename $i);
        dd=$(dirname $i)
        ff=$(basename $i)
        c1=`echo "$ff" | cut -d . -f 1`
        perl $data/bin/A9.1_Prepare_LEfSe.pl $i $dd/${ID%.*.*}.4lefse.txt $data/bin/$c1-color.txt 2;
done



source ~/miniconda2-py27/bin/activate 
source activate ~/miniconda2-py27/envs/lefse/
#Run LEfSe
for i in $(ls $data/Analysis/10_LEfSe/*/*.4lefse.txt);
do
	ID=$(basename $i);
	dd=$(dirname $i)
	lefse-format_input.py $i $dd/${ID%.*.*}.in -c 1 -u 2 -o 1000000;
	run_lefse.py $dd/${ID%.*.*}.in $dd/${ID%.*.*}.res;
	lefse-plot_res.py $dd/${ID%.*.*}.res --format pdf $dd/${ID%.*.*}.pdf;
	lefse-plot_features.py --format pdf --archive zip $dd/${ID%.*.*}.in $dd/${ID%.*.*}.res $dd/${ID%.*.*}.zip;
done

for i in $(ls $data/Analysis/10_LEfSe/*/*.4lefse.txt);
do
	ID=$(basename $i);
	dd=$(dirname $i)
	lefse-plot_cladogram.py $dd/${ID%.*.*}.res --format pdf $dd/${ID%.*.*}.clad.pdf --clade_sep 0.1;
	lefse-plot_cladogram.py $dd/${ID%.*.*}.res --format png $dd/${ID%.*.*}.clad.png --clade_sep 0.1;
done

