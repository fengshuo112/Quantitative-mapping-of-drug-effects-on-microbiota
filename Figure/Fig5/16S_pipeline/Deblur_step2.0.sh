data=$(dirname $PWD)

#[ ! -d $data/2_combined_fq ] && mkdir $data/2_combined_fq

time qiime cutadapt trim-paired \
--i-demultiplexed-sequences $data/1_import_data/demux.qza \
--p-front-f ACTCCTACGGGAGGCAGCAG \
--p-front-r GGACTACHVGGGTWTCTAAT  \
--o-trimmed-sequences $data/1_import_data/paired-end-demux.qza \
--verbose 
#&> primer_trimming.log
