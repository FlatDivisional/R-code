#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#引用包
library(limma)
setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\17.panGeneExp")    #设置工作目录
eRNA="AC003092.1"
taget="TFPI2"

#读取目录下的文件
files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
	#读取文件
	CancerType=gsub("symbol\\.|\\.txt","",i)
	rt=read.table(i,sep="\t",header=T,check.names=F)
	
	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	
	#得到样品的Type
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	Type=ifelse(group==0,"Tumor","Normal")
	outTab=rbind(outTab,cbind(t(data[c(eRNA,taget),]),Type,CancerType))
}
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="panGeneExp.txt",sep="\t",quote=F,row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
