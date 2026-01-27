# bash Deblur_step1_import_data.sh
echo step2
# bash Deblur_step2_Joined_fq.sh
echo step3
# bash Deblur_step3_filtered_fq.sh
echo step4
bash Deblur_step4_Picked_OTU.sh
echo step4.2
bash Deblur_step4.2_uchime.sh
echo step5
bash Deblur_step5.1_Modified_OTU.sh
bash Deblur_step5.2_Modified_Repseq.sh
echo step6
bash Deblur_step6_Phylogeny_Tree.sh
echo step7
bash Deblur_step7_Taxa_Greengene.sh
bash Deblur_step7_Taxa_Silva.sh
echo step8
bash Deblur_step8_Diversity.sh
bash Analysis_step0.sh
