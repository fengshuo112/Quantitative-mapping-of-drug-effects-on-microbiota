#Current working taxanomic level is Genus

Data=$(dirname $PWD)

##	MA/MAS/MC/MeH/MeS


for i in $(ls $Data/Analysis/15_RDA_and_Heatmap/*table.tsv)
do
	data=$(basename $i);
	d1=$(dirname $i)
	Out_dir=$Data/Analysis/15_RDA_and_Heatmap
	# up 和down 以及 生成物种信息
	arg3=("NC-PD" "NC-lhyy" "PD-lhyy" )  ##### A-B:B高为UP-1,B低为down-0,一般B为model
        Rscript $Data/bin/A11.3_Info_Factor.r $Out_dir ${arg2[*]} ${#arg1[*]} $Data/bin/Total-color.txt;
        #       factor:严格要求,同增同减
        # up 和down 以及 生成物种信息
        echo "A11.3"


	arg4=("NC-PD" "NC-lhyy" "PD-lhyy")  ####Heatmap 
	###### 不画grid
	Rscript $Data/bin/A11.4_Heatmap_no_grid.r $Out_dir ${arg4[*]}  $Data/bin/Total-color.txt;
	### 保存路径;分组信息;具体组
	##### 画 grid

	echo "A11.4"
	arg5=("PD"	"NC"	"lhyy") #Grid;第一位必须是基准即疾病;第二位必须是正常组，后面是给药组
	#Rscript $Data/bin/A11.4_Heatmap_grid.r $d1 ${arg4[*]}  $Data/bin/Total-color.txt ${arg5[*]} $p;
	### 保存路径;分组信息;具体组;grid信息




done



