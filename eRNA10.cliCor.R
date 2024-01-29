#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")

options(stringsAsFactors=F)
library(limma)
library(ggpubr)

setwd("E:\\BaiduNetdiskDownload\\eRNA����\\10.cliCor")     #�޸Ĺ���Ŀ¼
gene="TFPI2"           #eRNA���߰е�����

#��ȡ�����ļ������������ļ�����
rt=read.table("symbol.KIRC.txt",sep="\t",header=T,check.names=F) #���������ļ�������Ҫ�����о��������޸�
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ɾ��������Ʒ
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
data=rbind(data,gene=data[gene,])
exp=t(data[c("gene",gene),])

#��ȡ�ٴ������ļ�
cli=read.table("clinical.txt",sep="\t",header=T,check.names=F,row.names=1)
#cli[,"age"]=ifelse(cli[,"age"]=="unknow", "unknow", ifelse(cli[,"age"]>65,">65","<=65"))

#�ϲ�����
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
rt=cbind(exp,cli)

#�ٴ�����Է��������ͼ�ν��
for(clinical in colnames(rt[,3:ncol(rt)])){
	data=rt[c(gene,clinical)]
	colnames(data)=c("gene","clinical")
	data=data[(data[,"clinical"]!="unknow"),]
#���ñȽ���
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#����boxplot
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


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio