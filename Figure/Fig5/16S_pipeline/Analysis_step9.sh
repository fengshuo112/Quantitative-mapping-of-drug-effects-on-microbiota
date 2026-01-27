#!/bin/bash


Data=$(dirname $PWD)
PICRUST=/home/lijing/miniconda3/envs/picrust2/lib/python3.6/site-packages/picrust2/default_files

Rscript A10.0_part_kegg.r
[ ! -d $Data/Analysis/9_PICRUST2 ] && mkdir -p $Data/Analysis/9_PICRUST2

source activate picrust2
##无需创建文件夹。PICRUST会自己创建
picrust2_pipeline.py -s $Data/Analysis/0_input_data/rep-seqs.fasta \
        -i $Data/Analysis/0_input_data/Reotu_table.tsv \
        -o $Data/8_PICRUST_Result \
        -p 20
        
 # 生成kegg pathway 丰度表
pathway_pipeline.py -i$Data/8_PICRUST_Result/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \
    -o $Data/8_PICRUST_Result/KEGG_pathways_out --no_regroup \
    --map $PICRUST/pathway_mapfiles/KEGG_pathways_to_KO.tsv
    
# 添加功能描述
add_descriptions.py -i $Data/8_PICRUST_Result/KEGG_pathways_out/path_abun_unstrat.tsv.gz \
    --custom_map_table $PICRUST/description_mapfiles/KEGG_pathways_info.tsv.gz \
    -o $Data/8_PICRUST_Result/KEGG_pathways_out/path_abun_unstrat_descrip.tsv.gz

conda deactivate
 
 ##给KEGG通路加上二级标签
python3.8 $Data/bin/KEGG_Path_Change.py $Data/8_PICRUST_Result/KEGG_pathways_out/path_abun_unstrat_descrip.tsv.gz $Data/Analysis/9_PICRUST2/
 
Rscript $Data/bin/A3_Plot_PICRUST_box.r $Data/Analysis/9_PICRUST2/Total.catagrized.kegg.csv $Data/Analysis/9_PICRUST2/ $Data/bin/Total-color.txt

