#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")

#���ð�
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)


corFilter=0.4            #���ϵ�����˱�׼
pvalueFilter=0.001       #pֵ���˱�׼

setwd("C:\\Users\\lexb4\\Desktop\\eRNA\\09.corFilter")     #���ù���Ŀ¼
expFile="symbol.STAD.txt"            #���������ļ�,�ļ�������Ҫ�����о��������޸�

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

#��ȡ�������eRNA
eRNA=read.table("survResult.txt",sep="\t",header=T,check.names=F,row.names=1)

#��ȡTarget�ļ�
target=read.table("eRNAsymbol.txt",sep="\t",header=T,check.names=F)
target=target[,c(4,8)]
targetDup=target[!duplicated(target),]
targetDup=targetDup[which(targetDup[,1] %in% row.names(eRNA)),]
targetDup=targetDup[which(targetDup[,2] %in% row.names(data)),]

#����Լ���
outTab=data.frame()
for(n in 1:nrow(targetDup)){
	i=as.character(targetDup[n,1])
	j=as.character(targetDup[n,2])
	x=as.numeric(data[i,])
	y=as.numeric(data[j,])
	corT=cor.test(x,y,method="spearman")
	cor=corT$estimate
	pvalue=corT$p.value
	if((cor>corFilter) & (pvalue<pvalueFilter)){
		outTab=rbind(outTab,cbind(eRNA=i,KM=eRNA[i,],Target=j,cor=cor,corPval=pvalue))
		#�������������
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
			xlab(i)+ylab(j)+
			geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
		pdf(file=paste0("cor.",i,"_",j,".pdf"),width=5,height=5)
		print(p1)
		dev.off()
	}
}
#�������Խ��
write.table(file="corResult.xls",outTab,sep="\t",quote=F,row.names=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio