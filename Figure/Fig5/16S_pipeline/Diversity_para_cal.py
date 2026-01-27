import pandas as pd
import sys


def Calculate_filrate(dat,ra):
    Count_min = round(max(dat["count"])*ra,0)
    Satisifid_data = div.loc[dat["count"]>Count_min,]
    Satisifid_num = len(Satisifid_data)
    Total_num = len(dat)
    filter_rate = round(Satisifid_num/Total_num,3)
    return(filter_rate,Satisifid_data)

div = pd.read_csv(sys.argv[1],header=0,names=["sample","count"])
group = pd.read_table(sys.argv[2])


group_dict = {}
change_path = group.apply(lambda x:group_dict.setdefault(x["Group"],[]).append(x["sample-id"]),axis=1)

Difference_rate = round(min(div["count"])/max(div["count"]),3)
rate = Difference_rate if Difference_rate >0.1 else 0.1 #计算出极差比例，若大于0.1则保持不变，若小于0.1则使用0.1
filter_rate,data= Calculate_filrate(div,rate)
if filter_rate <0.5:
    rate = 0.05
    filter_rate,data = Calculate_filrate(div,rate)
else:
    pass

if filter_rate > 0.5:
    print("Your data is Ok!")
    output1 = "Count:"+str(round(max(div["count"])*rate,0))+"\n"
    output2 = "Count_rate:"+str(rate)+"\n"
    output3 = "Filter_rate:"+str(filter_rate)+"\n"
    output = list(output1+output2+output3)
    with open(sys.argv[3],"w+") as f:
        f.writelines(output)
else:
    print("Your data is uneven distibution!")
    output1 = "Count:"+str(round(max(div["count"])*rate,0))+"\n"
    output2 = "Count_rate"+":"+str(rate)+"\n" 
    output3 = "Filter_rate"+":"+str(filter_rate)+"\n"
    output = list(output1+output2+output3)
    with open(sys.argv[3],"w+") as f:
        f.writelines(output)