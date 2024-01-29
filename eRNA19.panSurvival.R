#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

gene="TFPI2"       #��������
pFilter=0.05        #km����pvalue��������
setwd("E:\\BaiduNetdiskDownload\\eRNA����\\19.panSurvival")      #���ù���Ŀ¼
rt=read.table("panExpTime.txt",header=T,sep="\t",check.names=F,row.names=1)    #��ȡ�����ļ�
rt$futime=rt$futime/365

#���������ͽ���ѭ��
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
		#������������
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
				    width = 6,             #ͼƬ�Ŀ���
				    height =5)             #ͼƬ�ĸ߶�
		print(surPlot)
		dev.off()
	}
}
#����������ͺ�pֵ�����ļ�
write.table(outTab,file="panSurvResult.xls",sep="\t",row.names=F,quote=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio