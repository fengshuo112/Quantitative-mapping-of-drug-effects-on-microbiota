#Current working taxanomic level is Genus

Data=$(dirname $PWD)

##	MA/MAS/MC/MeH/MeS
[ ! -d $Data/Analysis/15_RDA_and_Heatmap ] && mkdir -p $Data/Analysis/15_RDA_and_Heatmap

cp $Data/Analysis/0_input_data/Reotu_table.tsv $Data/Analysis/15_RDA_and_Heatmap
cp $Data/Analysis/0_input_data/taxonomy.tsv $Data/Analysis/15_RDA_and_Heatmap


for i in $(ls $Data/Analysis/15_RDA_and_Heatmap/*table.tsv)
do
	data=$(basename $i);
	d1=$(dirname $i)
	Out_dir=$Data/Analysis/15_RDA_and_Heatmap
	arg1=("B_ctl-M_ctl" "B_ctl-M_D" "B_ctl-M_E" "M_ctl-M_D" "M_ctl-M_E" "M_D-M_E" )	####分组信息
	Rscript $Data/bin/A11.1_RDA.r $data $d1/taxonomy.tsv $Out_dir ${arg1[*]}  $Data/bin/Total-color.txt;
	echo "A11.1"
	##参数为读入文件;保存路径;分组信息
	###### 主要生成p值
	p=0.05

	arg2=("B_ctl-M_ctl" "B_ctl-M_D" "B_ctl-M_E" "M_ctl-M_D" "M_ctl-M_E" "M_D-M_E")  #####
	Rscript $Data/bin/A11.2_Info_P.r $Out_dir ${arg2[*]} $p  ${#arg1[*]} $Data/bin/Total-color.txt;	
	#参数为保存路径;分组信息;p值;分组数;原始的分组信息
	####### 挑选p值,全部小于P 
	echo "A11.2"

	arg3=("B_ctl-M_ctl" "B_ctl-M_D" "B_ctl-M_E" "M_ctl-M_D" "M_ctl-M_E" "M_D-M_E" )  ##### A-B:B高为UP-1,B低为down-0,一般B为model
	Rscript $Data/bin/A11.3_Info_Factor.r $Out_dir ${arg2[*]} ${#arg1[*]} $Data/bin/Total-color.txt;	
	#	factor:严格要求,同增同减
	# up 和down 以及 生成物种信息
	echo "A11.3"


	arg4=("B_ctl" "M_ctl" "M_D" "M_E" )  ####Heatmap 
	###### 不画grid
	Rscript $Data/bin/A11.4_Heatmap_no_grid.r $Out_dir ${arg4[*]}  $Data/bin/Total-color.txt;
	### 保存路径;分组信息;具体组
	##### 画 grid

	echo "A11.4"
	arg5=("PD"	"NC"	) #Grid;第一位必须是基准即疾病;第二位必须是正常组，后面是给药组
#	Rscript $Data/bin/A11.4_Heatmap_grid.r $d1 ${arg4[*]}  $Data/bin/Total-color.txt ${arg5[*]} $p;

	### 保存路径;分组信息;具体组;grid信息




done



