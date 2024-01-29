#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)     #���ð�
setwd("E:\\BaiduNetdiskDownload\\eRNA����\\11.cliStat")     #�޸Ĺ���Ŀ¼

#��ȡ�����ļ������������ļ�����
rt=read.table("symbol.KIRC.txt",sep="\t",header=T,check.names=F)  #���������ļ�������Ҫ�����о��������޸�
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

#��ȡ�ٴ������ļ�
cli=read.table("clinical.txt",sep="\t",header=T,check.names=F,row.names=1)

#�ٴ����ݺͱ�������ȡ����
samSample=intersect(colnames(data),row.names(cli))
cli=cli[samSample,]

#���ٴ����ݽ���ͳ��
cliStatOut=data.frame()
for(i in 1:ncol(cli)){
	nameStat=colnames(cli)[i]
	tableStat=table(cli[,i])
	tableStatSum=cbind(Total=sum(tableStat),tableStat)
	tableStatRatio=prop.table(tableStatSum,2)
	tableStatRatio=round(tableStatRatio*100,2)
	tableStatSum[,2]=paste(tableStatSum[,2],"(",tableStatRatio[,2],"%)",sep="")
	tableStatOut=cbind(Covariates=nameStat,Type=row.names(tableStatSum),tableStatSum)
	cliStatOut=rbind(cliStatOut,tableStatOut)
}
#����ٴ�ͳ�ƽ��
cliStatOut=cliStatOut[,-3]
write.table(cliStatOut,file="cliStat.xls",sep="\t",quote=F,row.names=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio