#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#���ð�
library(limma)

expFile="panGeneExp.txt"                                        #���������ļ�
setwd("E:\\BaiduNetdiskDownload\\eRNA����\\18.panMergeTime")       #���ù���Ŀ¼

#��ȡĿ¼�����������ļ�
files=dir()
files=grep(".survival.tsv$",files,value=T)

#��ȡÿ�����������������ļ�
surTime=data.frame()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime,rt)
}
colnames(surTime)=c("futime","fustat")

#��ȡ�����ļ������������ļ�����
exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

#���ݺϲ���������
sameSample=intersect(row.names(surTime),row.names(exp))
surTime=surTime[sameSample,]
exp=exp[sameSample,]
surData=cbind(surTime,exp)
surData=cbind(id=row.names(surData),surData)
write.table(surData,file="panExpTime.txt",sep="\t",quote=F,row.names=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio