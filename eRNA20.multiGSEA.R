#install.packages("ggplot2")

library(plyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("E:\\BaiduNetdiskDownload\\30.multiGSEA")          #设置工作目录
files=grep(".xls",dir(),value=T)                                     #获取目录下的所有xls文件
data = lapply(files,read.delim)                                      #读取每个文件
names(data) = files

dataSet = ldply(data, data.frame)
dataSet$pathway = gsub(".xls","",dataSet$.id)                            #将文件后缀删掉

gseaCol=c("#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
pGsea=ggplot(dataSet,aes(x=RANK.IN.GENE.LIST,y=RUNNING.ES,colour=pathway,group=pathway))+
  geom_line(size = 1.5) + scale_color_manual(values = gseaCol[1:nrow(dataSet)]) +   
  labs(x = "", y = "Enrichment Score", title = "") + scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0),limits = c(min(dataSet$RUNNING.ES - 0.02), max(dataSet$RUNNING.ES + 0.02))) +   
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank()) + 
  geom_hline(yintercept = 0) + 
  guides(colour = guide_legend(title = NULL)) + theme(legend.background = element_blank()) + theme(legend.key = element_blank())+theme(legend.key.size=unit(0.5,'cm'))
pGene=ggplot(dataSet,aes(RANK.IN.GENE.LIST,pathway,colour=pathway))+geom_tile()+
  scale_color_manual(values = gseaCol[1:nrow(dataSet)]) + 
  labs(x = "high expression<----------->low expression", y = "", title = "") + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  
  theme_bw() + theme(panel.grid = element_blank()) + theme(panel.border = element_blank()) + theme(axis.line = element_line(colour = "black"))+
  theme(axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+ guides(color=FALSE)

gGsea = ggplot_gtable(ggplot_build(pGsea))
gGene = ggplot_gtable(ggplot_build(pGene))
maxWidth = grid::unit.pmax(gGsea$widths, gGene$widths)
gGsea$widths = as.list(maxWidth)
gGene$widths = as.list(maxWidth)
dev.off()

#将图形可视化，保存在"multipleGSEA.pdf"
pdf('multipleGSEA.pdf',      #输出图片的文件
     width=10,               #设置输出图片高度
     height=6)               #设置输出图片高度
par(mar=c(5,5,2,5))
grid.arrange(arrangeGrob(gGsea,gGene,nrow=2,heights=c(.8,.3)))
dev.off()


######生信自学网: http://study.163.com/u/biowolf
######生信自学网: https://shop119322454.taobao.com
######生信自学网: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
