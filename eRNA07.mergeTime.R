#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#引用包
library(limma)

expFile="eRNAexp.txt"                                       #表达数据文件
surFile="TCGA-KIRC.survival.tsv"                            #生存数据文件,文件名称需要根据研究的肿瘤修改
setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\07.mergeTime")      #设置工作目录

#读取生存数据
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
surTime=surTime[,c(3,1)]
colnames(surTime)=c("futime","fustat")

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
