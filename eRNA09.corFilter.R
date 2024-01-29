#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")

#引用包
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)


corFilter=0.4            #相关系数过滤标准
pvalueFilter=0.001       #p值过滤标准

setwd("C:\\Users\\lexb4\\Desktop\\eRNA\\09.corFilter")     #设置工作目录
expFile="symbol.STAD.txt"            #表达数据文件,文件名称需要根据研究的肿瘤修改

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

#读取生存相关eRNA
eRNA=read.table("survResult.txt",sep="\t",header=T,check.names=F,row.names=1)

#读取Target文件
target=read.table("eRNAsymbol.txt",sep="\t",header=T,check.names=F)
target=target[,c(4,8)]
targetDup=target[!duplicated(target),]
targetDup=targetDup[which(targetDup[,1] %in% row.names(eRNA)),]
targetDup=targetDup[which(targetDup[,2] %in% row.names(data)),]

#相关性检验
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
		#绘制相关性曲线
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
#输出相关性结果
write.table(file="corResult.xls",outTab,sep="\t",quote=F,row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
