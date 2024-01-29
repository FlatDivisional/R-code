#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)
setwd("E:\\BaiduNetdiskDownload\\eRNA����\\12.geneCor")     #���ù���Ŀ¼
expFile="symbol.KIRC.txt"                                #�����ļ�,�ļ�������Ҫ�����о��������޸�
gene="AC003092.1"                                            #�����lncRNA����
corFilter=0.4                                            #���ϵ������ֵ
pFilter=0.001                                            #ͳ��ѧpֵ����ֵ

#��ȡ�����ļ������������ļ�����
rt=read.table(expFile,sep="\t",header=T,check.names=F)
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

    #�������Խ��
	if((cor>corFilter) & (pvalue<pFilter)){
		outTab=rbind(outTab,cbind(eRNA=gene1,gene=gene2,cor,pvalue))
	}
}
#�������Խ��
outputFile=paste(gene1,".cor.xls",sep="")
write.table(file=outputFile,outTab,sep="\t",quote=F,row.names=F)
outputFile=paste(gene1,".cor.txt",sep="")
write.table(file=outputFile,outTab,sep="\t",quote=F,row.names=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio