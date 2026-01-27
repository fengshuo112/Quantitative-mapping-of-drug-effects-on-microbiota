arg1 <- commandArgs(T)
if(length(arg1) < 0){
        cat("Argument: Data_File Group_File Out_Dir Group_column\n")
        quit('no')
}
print(arg1)
# arg1=c(  "/home/wmj/视频/课题16s/tt/QC_Raw/Result_Table_Info_Raw.txt"  )


################ Raw
data=read.table(arg1[1],sep='\t',header=T,row.names=1)
d2=as.numeric(data[,3])
mm=ceiling(mean(d2))
print(paste0('AvgLen_Raw_Pair : ',mm))

