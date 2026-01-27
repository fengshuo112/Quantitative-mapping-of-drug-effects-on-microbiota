from numpy import float64
import pandas as pd
from Bio.KEGG import REST
import gzip
import sys

print(sys.argv)

def Usage():
    print("Python KEGG_Path_Change.py Input_File Output_Dir")

if len(sys.argv) <2:
    Usage()
else:
    with gzip.open(sys.argv[1],"rb") as f:
        data = f.readlines()
    data_list = list(map(lambda x:x.decode().strip().split("\t"),data))
    df = pd.DataFrame(data_list[1:],columns=data_list[0],dtype=float64)

    pwid = df['pathway']
    my_dict={}
    list_nofind= []
    for i,id in enumerate(pwid):
        try:
            pathway = REST.kegg_get(id).read()
        except:
            list_nofind.append(id)
        #print(i,id)
        for i in pathway.rstrip().split('\n'):
            feature = str(i.split('\t'))
            if 'CLASS 'in feature:
                q = str(i.split(';')[0])
                name1 = str(q.strip('CLASS').strip())
    #            name2 = str(i.split(';')[1])
    #            value = name1 +'/'+ name2
                my_dict[id] = name1
    #             print(my_dict)

    df["abun"] = df.apply(lambda x:round(x[2:].mean(),0),axis=1)
    data_all = df.loc[:,"abun"].sum()
    df["relati.abun"] = df["abun"].apply(lambda x:round(x/data_all,3))
    df['name'] = df['pathway'].apply(lambda x : my_dict.get(x,"nofind"))

    df.to_csv(sys.argv[2]+"Total.catagrized.kegg.csv",index=None)

