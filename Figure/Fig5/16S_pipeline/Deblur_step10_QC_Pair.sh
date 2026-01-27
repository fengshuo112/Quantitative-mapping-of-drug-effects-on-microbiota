### 双端序列合并
Data=$(dirname $PWD)
cd $Data;


###### raw
[ ! -d QC_Raw ] && mkdir QC_Raw;
for i in $(ls $Data/0_raw_data/*.fastq); 
do
	fastqc -o $Data/QC_Raw -t 40 $i --extract ## 一般几个样品就用几个线程
done
perl $Data/bin/A5.1_Extract_QC_Raw.pl $Data'/0_raw_data/*.fastq' $Data'/QC_Raw/*_fastqc/fastqc_data.txt' $Data/QC_Raw
rm -r $Data/QC_Raw/*.zip

Rscript $Data/bin/A5_QC_Avelen.r $Data"/QC_Raw/Result_Table_Info_Raw.txt" 
####### 参数说明
#	-o，结果输出路径；
#	--extract，默认情况下，会将所有结果文件打包后生成一个压缩文件，该参数存在时将解压缩；
#	-t，程序运行时的线程数；
#	-c，指定一个文件，该文件内容为“name[tab]sequence”样式，记录可能的污染序列， FastQC会根据该文件中的序列信息评估测序数据的污染程度，不指定时不评估；
#	-a，指定一个文件，该文件内容为“name[tab]sequence”样式，记录测序接头序列，FastQC会根据该文件中的序列信息评估测序接头序列的残留情况，不指定时将自动识别通用接头序列；
#	-k，指定k-mer统计时的k-mer长度，取值范围2-10，默认7；
#	-q，默认情况下，程序会实时报告运行的状况，该参数存在时仅报告错误信息。


######## merge
[ ! -d QC_Merge ] && mkdir QC_Merge;
for i in $(ls $Data/2_combined_fq/*.fastq); 
do
	fastqc -o $Data/QC_Merge -t 4 $i --extract ## 一般几个样品就用几个线程
done
perl $Data/bin/A5.2_Extract_QC_Merge.pl $Data'/2_combined_fq/*.fastq' $Data'/QC_Merge/*_fastqc/fastqc_data.txt' $Data/QC_Merge
rm -r $Data/QC_Merge/*.zip



[ ! -d QC_End ] && mkdir QC_End;
perl $Data/bin/A5.3_Extract_QC_End.pl   $Data/3_filter_fq/"*.fa" $Data/QC_End 


Rscript $Data/bin/A5.4_QC_Plot_Pair.r $Data"/QC_Raw/*LengthCount.txt" $Data"/QC_Raw/Result_Table_Info_Raw.txt" $Data"/QC_Merge/Result_Table_Info_Merge.txt" $Data/QC_End/Result_Length_End.txt $Data/QC_End/Result_Table_Info_End.txt  $Data/QC_End 
 

