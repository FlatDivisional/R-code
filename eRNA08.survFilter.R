#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\08.survFilter")                      #工作目录（需修改）
pFilter=0.01                                                                 #显著性过滤条件
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
rt$futime=rt$futime/365                                                      #除以365，把单位改成年

#对eRNA进行循环，做生存分析
outTab=data.frame()
for(gene in colnames(rt[,3:ncol(rt)])){
	#删除表达波动小的eRNA
	if(sd(rt[,gene])<0.1){
	   next}
	#根据eRNA中位值对样品分组
	a=ifelse(rt[,gene]<=median(rt[,gene]),"Low","High")
	#高低组生存差异
	diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	#计算每个时间节点病人数目
	fit=survfit(Surv(futime, fustat) ~ a, data = rt)
    #对满足条件的eRNA绘制生存曲线
	if(pValue<pFilter){
	    outTab=rbind(outTab,cbind(gene=gene,KM=pValue) )
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=sprintf("%.03f",pValue)
			pValue=paste0("p=",pValue)
		}
		titleName=paste0(gene," level")
		surPlot=ggsurvplot(fit, 
							data=rt,
							conf.int=TRUE,
							pval=pValue,
							pval.size=6,
							risk.table=T,
							legend.labs=c("high","low"),
							legend.title=titleName,
							xlab="Time(years)",
							break.time.by = 1,
							risk.table.title="",
							palette=c("red", "blue"),
							risk.table.height=.25)          
		pdf(file=paste0("sur.",gene,".pdf"), width = 6.5, height = 5.5,onefile = FALSE)
		print(surPlot)
		dev.off()
	}
}
#输出基因和p值表格文件
write.table(outTab,file="survResult.txt",sep="\t",row.names=F,quote=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
