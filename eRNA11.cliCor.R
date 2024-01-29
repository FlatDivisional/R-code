#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)     #引用包
setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\11.cliStat")     #修改工作目录

#读取表达文件，并对输入文件整理
rt=read.table("symbol.KIRC.txt",sep="\t",header=T,check.names=F)  #表达数据文件名称需要根据研究的肿瘤修改
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

#读取临床数据文件
cli=read.table("clinical.txt",sep="\t",header=T,check.names=F,row.names=1)

#临床数据和表达数据取交集
samSample=intersect(colnames(data),row.names(cli))
cli=cli[samSample,]

#对临床数据进行统计
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
#输出临床统计结果
cliStatOut=cliStatOut[,-3]
write.table(cliStatOut,file="cliStat.xls",sep="\t",quote=F,row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
