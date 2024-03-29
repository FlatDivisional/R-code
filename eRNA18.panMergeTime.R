#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#引用包
library(limma)

expFile="panGeneExp.txt"                                        #表达数据文件
setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\18.panMergeTime")       #设置工作目录

#获取目录下所有生存文件
files=dir()
files=grep(".survival.tsv$",files,value=T)

#读取每个肿瘤的生存数据文件
surTime=data.frame()
for(i in 1:length(files)){
    inputFile=files[i]
    rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
    rt=rt[,c(3,1)]
    surTime=rbind(surTime,rt)
}
colnames(surTime)=c("futime","fustat")

#读取表达文件，并对输入文件整理
exp=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
exp=exp[(exp[,"Type"]=="Tumor"),]

#数据合并并输出结果
sameSample=intersect(row.names(surTime),row.names(exp))
surTime=surTime[sameSample,]
exp=exp[sameSample,]
surData=cbind(surTime,exp)
surData=cbind(id=row.names(surData),surData)
write.table(surData,file="panExpTime.txt",sep="\t",quote=F,row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
