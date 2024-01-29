#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

gene="TFPI2"       #基因名称
pFilter=0.05        #km方法pvalue过滤条件
setwd("E:\\BaiduNetdiskDownload\\eRNA资料\\19.panSurvival")      #设置工作目录
rt=read.table("panExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)    #读取输入文件
rt$futime=rt$futime/365

#对肿瘤类型进行循环
outTab=data.frame()
for(i in levels(factor(rt[,"CancerType"]))){
	rt1=rt[(rt[,"CancerType"]==i),]
	group=ifelse(rt1[,gene]>median(rt1[,gene]),"high","low")
	diff=survdiff(Surv(futime, fustat) ~group,data = rt1)
	pValue=1-pchisq(diff$chisq,df=1)
	outTab=rbind(outTab,cbind(gene=gene,CancerType=i,KM=pValue) )
	if(pValue<pFilter){
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
		titleName=paste0(gene," level")
		#绘制生存曲线
		surPlot=ggsurvplot(fit, 
					data=rt1,
					title=paste0("Cancer: ",i),
					conf.int=TRUE,
					pval=pValue,
					pval.size=6,
					risk.table=T,
					legend.labs=c("high","low"),
					legend.title=titleName,
					xlab="Time(years)",
					break.time.by = 1,
					risk.table.title="",
					palette=c("red", "blue"))
		pdf(file=paste0("survival.",i,".pdf"),onefile = FALSE,
				    width = 6,             #图片的宽度
				    height =5)             #图片的高度
		print(surPlot)
		dev.off()
	}
}
#输出肿瘤类型和p值表格文件
write.table(outTab,file="panSurvResult.xls",sep="\t",row.names=F,quote=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱：seqbio@foxmail.com
######光俊老师微信: seqBio
