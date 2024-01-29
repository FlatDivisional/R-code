#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\12.geneCor")     #设置工作目录
expFile="symbol.KIRC.txt"                                #输入文件,文件名称需要根据研究的肿瘤修改
gene="AC003092.1"                                            #基因或lncRNA名字
corFilter=0.4                                            #相关系数过滤值
pFilter=0.001                                            #统计学p值过滤值

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

x=as.numeric(data[gene,])
gene1=unlist(strsplit(gene,"\\|",))[1]

outTab=data.frame()
for(j in rownames(data)){
    y=as.numeric(data[j,])
    if(sd(y)<0.1){
    	next
    }
    if(gene==j){
    	next
    }
	gene2=unlist(strsplit(j,"\\|",))[1]
	corT=cor.test(x,y,method="spearman")
	cor=corT$estimate
	cor=round(cor,3)
	pvalue=corT$p.value

    #输出相关性结果
	if((cor>corFilter) & (pvalue<pFilter)){
		outTab=rbind(outTab,cbind(eRNA=gene1,gene=gene2,cor,pvalue))
	}
}
#输出相关性结果
outputFile=paste(gene1,".cor.xls",sep="")
write.table(file=outputFile,outTab,sep="\t",quote=F,row.names=F)
outputFile=paste(gene1,".cor.txt",sep="")
write.table(file=outputFile,outTab,sep="\t",quote=F,row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
