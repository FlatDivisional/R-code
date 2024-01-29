#install.packages("survival")
#install.packages("survminer")

library(survival)
library(survminer)

setwd("E:\\BaiduNetdiskDownload\\eRNA����\\08.survFilter")                      #����Ŀ¼�����޸ģ�
pFilter=0.01                                                                 #�����Թ�������
rt=read.table("expTime.txt",header=T,sep="\t",check.names=F,row.names=1)     #��ȡ�����ļ�
rt$futime=rt$futime/365                                                      #����365���ѵ�λ�ĳ���

#��eRNA����ѭ�������������
outTab=data.frame()
for(gene in colnames(rt[,3:ncol(rt)])){
	#ɾ�����ﲨ��С��eRNA
	if(sd(rt[,gene])<0.1){
	   next}
	#����eRNA��λֵ����Ʒ����
	a=ifelse(rt[,gene]<=median(rt[,gene]),"Low","High")
	#�ߵ����������
	diff=survdiff(Surv(futime, fustat) ~a,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	#����ÿ��ʱ��ڵ㲡����Ŀ
	fit=survfit(Surv(futime, fustat) ~ a, data = rt)
    #������������eRNA������������
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
#��������pֵ�����ļ�
write.table(outTab,file="survResult.txt",sep="\t",row.names=F,quote=F)


######������ѧ��: https://www.biowolf.cn/
######�γ�����1: https://shop119322454.taobao.com
######�γ�����2: https://ke.biowolf.cn
######�γ�����3: https://ke.biowolf.cn/mobile
######�⿡��ʦ���䣺seqbio@foxmail.com
######�⿡��ʦ΢��: seqBio