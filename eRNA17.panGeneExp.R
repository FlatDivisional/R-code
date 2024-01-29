#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#���ð�
library(limma)
setwd("E:\\BaiduNetdiskDownload\\eRNA����\\17.panGeneExp")    #���ù���Ŀ¼
eRNA="AC003092.1"
taget="TFPI2"

#��ȡĿ¼�µ��ļ�
files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
	#��ȡ�ļ�
	CancerType=gsub("symbol\\.|\\.txt","",i)
	rt=read.table(i,sep="\t",header=T,check.names=F)
	
	#���һ������ռ�˶��У�ȡ��ֵ
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	
	#�õ���Ʒ��Type
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	Type=ifelse(group==0,"Tumor","Normal")
	outTab=rbind(outTab,cbind(t(data[c(eRNA,taget),]),Type,CancerType))
}
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="panGeneExp.txt",sep="\t",quote=F,row.names=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio