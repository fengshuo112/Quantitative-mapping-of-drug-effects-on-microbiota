#!/bin/bash

Data=$(dirname $PWD)

#ls $Data/0_raw_data/*fastq.gz|perl -e 'while(<>){chomp;@a=split/\//;@b=split/\_/,$a[-1];$inf=<>;print"$b[0]\t$_\t$inf"}'|awk 'OFS="\t"{if($0~/x/) print $0,"x","#D4D4D4";else print $0,"y","#1874CD"}'|sed '1i sample-id\tforward-absolute-filepath\treverse-absolute-filepath\tGroup\tColour' > sample-metadata.txt

#ls $Data/0_raw_data/*fq.gz|perl -e 'while(<>){chomp;@a=split/\//;@b=split/\./,$a[-1];print"$b[0]\t$_\n"}'|awk 'OFS="\t"{if($0~/1/) print $0,"Z1","#D4D4D4";else print $0,"Z2","#1874CD"}'|sed '1i sample-id\tabsolute-filepath\tGroup\tColour' > sample-metadata.txt


[ ! -d $Data/1_import_data ] && mkdir -p $Data/1_import_data
[ ! -d $Data/1_import_data/data ] && mkdir -p $Data/1_import_data/data

#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
#	--input-path $Data/bin/sample-metadata.txt \
#	--output-path $Data/1_import_data/demux.qza \
#	--input-format PairedEndFastqManifestPhred33V2

qiime tools import --type 'SampleData[SequencesWithQuality]' \
	--input-path $Data/bin/sample-metadata.txt \
	--output-path $Data/1_import_data/demux.qza \
	--input-format SingleEndFastqManifestPhred33V2
