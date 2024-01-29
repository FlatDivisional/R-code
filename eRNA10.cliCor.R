#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

options(stringsAsFactors=F)
library(limma)
library(ggpubr)

setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\10.cliCor")     #修改工作目录
gene="TFPI2"           #eRNA或者靶点名称

#读取表达文件，并对输入文件整理
rt=read.table("symbol.KIRC.txt",sep="\t",header=T,check.names=F) #表达数据文件名称需要根据研究的肿瘤修改
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
data=rbind(data,gene=data[gene,])
exp=t(data[c("gene",gene),])

#读取临床数据文件
cli=read.table("clinical.txt",sep="\t",header=T,check.names=F,row.names=1)
#cli[,"age"]=ifelse(cli[,"age"]=="unknow", "unknow", ifelse(cli[,"age"]>65,">65","<=65"))

#合并数据
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
rt=cbind(exp,cli)

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,3:ncol(rt)])){
	data=rt[c(gene,clinical)]
	colnames(data)=c("gene","clinical")
	data=data[(data[,"clinical"]!="unknow"),]
#设置比较组
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#绘制boxplot
	boxplot=ggboxplot(data, x="clinical", y="gene", color="clinical",
	          xlab=clinical,
	          ylab=paste(gene,"expression"),
	          legend.title=clinical,
	          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	pdf(file=paste0(clinical,".pdf"),width=5.5,height=5)
	print(boxplot)
	dev.off()
}


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
