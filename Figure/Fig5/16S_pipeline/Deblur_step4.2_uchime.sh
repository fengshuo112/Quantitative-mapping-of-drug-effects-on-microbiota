Data=$(dirname $PWD)

time qiime vsearch uchime-denovo \
  --i-table $Data/4_Picked_OTUs/table-99.qza \
  --i-sequences $Data/4_Picked_OTUs/rep-seqs-99.qza \
  --output-dir $Data/4_Picked_OTUs/uchime-dn-out

time qiime metadata tabulate \
  --m-input-file $Data/4_Picked_OTUs/uchime-dn-out/stats.qza \
  --o-visualization $Data/4_Picked_OTUs/uchime-dn-out/stats.qzv

time qiime feature-table filter-features \
  --i-table $Data/4_Picked_OTUs/table-99.qza \
  --m-metadata-file $Data/4_Picked_OTUs/uchime-dn-out/nonchimeras.qza \
  --o-filtered-table $Data/4_Picked_OTUs/uchime-dn-out/table-nonchimeric-wo-borderline.qza
time qiime feature-table filter-seqs \
  --i-data $Data/4_Picked_OTUs/rep-seqs-99.qza \
  --m-metadata-file $Data/4_Picked_OTUs/uchime-dn-out/nonchimeras.qza \
  --o-filtered-data $Data/4_Picked_OTUs/uchime-dn-out/rep-seqs-nonchimeric-wo-borderline.qza
time qiime feature-table summarize \
  --i-table $Data/4_Picked_OTUs/uchime-dn-out/table-nonchimeric-wo-borderline.qza \
  --o-visualization $Data/4_Picked_OTUs/uchime-dn-out/table-nonchimeric-wo-borderline.qzv
