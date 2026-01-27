#Current working taxanomic level is Genus

#Data=/home/wmj/Metagenome-2/16s/ZSL-1
Data=$(dirname $PWD)
cd $Data;

[ ! -d $Data/Analysis/17_kegg_L3 ] && mkdir $Data/Analysis/17_kegg_L3;


for i in $(ls $Data/Analysis/10_LEfSe/*/*kegg.L3.res);
do
	ff=$(basename $i)
	dd=$(dirname $i)
	d1=${dd/10_LEfSe/17_kegg_L3};
	[ ! -d $d1 ] && mkdir $d1;
	cp $i $d1

	File=$(basename $i)
	c1=`echo "$File" | cut -d . -f 1`
	Rscript $Data/bin/A13_LDA_plot.r $i $d1 $Data/bin/$c1-color.txt

done

for i in $(ls $Data/Analysis/10_LEfSe/*/*Genus.res);
do
        ff=$(basename $i)
        dd=$(dirname $i)
        d1=${dd/10_LEfSe/17_kegg_L3};
        [ ! -d $d1 ] && mkdir $d1;
        cp $i $d1

        File=$(basename $i)
        c1=`echo "$File" | cut -d . -f 1`
        Rscript $Data/bin/A13_LDA_plot.r $i $d1 $Data/bin/$c1-color.txt

done
