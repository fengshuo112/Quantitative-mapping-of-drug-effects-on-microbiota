### 双端序列合并
Data=/home/Project/16S/20200225_WXD
cd $Data;


###### raw
[ ! -d QC_Raw ] && mkdir QC_Raw;
for i in $(ls $Data/1_combine_fq/*.fq); 
do
	fastqc -o $Data/QC_Raw -t 4 $i --extract ## 一般几个样品就用几个线程
done
perl $Data/bin/A5.1_Extract_QC_Raw.pl $Data'/1_combine_fq/*.fq' $Data'/QC_Raw/*_fastqc/fastqc_data.txt' $Data/QC_Raw
rm -r $Data/QC_Raw/*.zip

Rscript $Data/bin/A5_QC_Avelen.r $Data"/QC_Raw/Result_Table_Info_Raw.txt" 
 


[ ! -d QC_End ] && mkdir QC_End;
perl $Data/bin/A5.3_Extract_QC_End.pl   $Data'/3_filtered_fa/*.fa' $Data/QC_End 


Rscript $Data/bin/A5.4_QC_Plot_Single.r $Data"/QC_Raw/*LengthCount.txt" $Data"/QC_Raw/Result_Table_Info_Raw.txt"  $Data/QC_End/Result_Length_End.txt $Data/QC_End/Table.txt  $Data/QC_End  
cp $Data"/QC_Raw/Result_Table_Info_Raw.txt" $Data"/QC_End/Table.txt" 
