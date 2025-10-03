library(ggplot2)
library(eoffice)
#library(extrafont)
library(export)
library(rtracklayer)
library(GenomicFeatures)
library(movAPA)
library(wesanderson)
library(dplyr)
library(tidyr)
library(patchwork)
library(aplot)
library(RColorBrewer)
library(pROC)
library(UpSetR)
library(psych)
library(ggpubr)
library(cowplot) 
library(SummarizedExperiment)
library(Seurat)
library(cluster)


tool.col=c("scAPAtrap"="#064D4B","Sierra"="#fbd26a","SCAPE"="#ABC8E5","scAPA"="#6666CC","polyApipe"="#e8b9c5","MAAPER"="#64d9b8",
           "SCAPTURE"="#b5182b","Infernape"="#775360","scUTRquant"="#7fb80e","scraps"="#FF6B00")

p_theme<- theme_bw()+theme(
  text=element_text(family="Arial" ),
  axis.title.x =element_text(color="black",size=10,family="Arial" ) ,
  axis.title.y =element_text(color="black",size=10,family="Arial" ) ,
  axis.text.x =element_text(color="black",size=8,family="Arial" ) ,
  axis.text.y =element_text(color="black",size=9,family="Arial" ) ,
  legend.text =element_text(color="black",size=8,family="Arial"),
  legend.title=element_text(color="black",size=8,family= "Arial" ),
  legend.background = element_blank(),
  panel.border = element_blank(),panel.grid.major = element_blank(),
  panel.grid.minor =element_blank(),
  axis.line=element_line(colour = "black",size=0.4))

###1.Comparison of time and memory consumption of different methods----------------
res=readxl::read_excel("G:/scAPA/summary/temp/tool.ues.resource.xlsx")
res <- res %>%
  mutate(
    `Run Time (Min)` = as.numeric(`Run Time (Min)`),
    `Max memory (MB)` = as.numeric(`Max memory (MB)`)
  )


res$Method=factor(res$Method,
                  levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
#`Run Time (Min)` `Max memory (MB)`


###(1.1)run time--------------------------------------------------
custom_breaks=unique(res$`Cell Number`)

p=ggplot(res, aes(x = `Cell Number`,y =`Run Time (Min)`/60 , color =Method,fill= Method,group = Method)) +  #group=group：确保 geom_line() 按照 group 分组连接点
  #geom_line(linetype="dashed",size=1) +
  geom_line() +
  geom_point(size=1.8,shape=23,
             stroke = 0.5,           
             alpha=0.9) +  
  scale_color_manual(values=tool.col)+
  scale_fill_manual(values=tool.col)+
  labs(x = "Cell Number", y = "Run Time (Hour)")+
  p_theme+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1,size=7),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    # 隐藏断点区域的 X 轴文本
    axis.text.x.top = element_blank(),     
    axis.ticks.x.top = element_blank(),   
    axis.line.x.top = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.text = element_text(size = 7)
   
  )+
  #guides(fill="none")+
  scale_x_continuous(breaks = custom_breaks,  
                     labels = custom_breaks, 
                     limits = range(0,max(custom_breaks)))+ 
  guides(color=guide_legend(title=''),fill='none')+ 
  ggbreak::scale_y_break(
    #breaks =c(8000, 11000),  # 截断区间
    breaks =c(120, 150),
    space = 0.1,      
    ticklabels=c(150,530), 
    scales = 0.25
   
  ) 
  
  #theme(plot.margin = margin(t = 80, unit = "pt")) 
  # 第二个截断区间
  #ggbreak::scale_y_break(
  #  breaks = c(11500, 30000), 
  #  space = 0.1,
  #  ticklabels=c(30000,31600), 
  #  scales = "fixed")
  


p1=ggplotify::as.ggplot(p) 


topptx(p1, filename = paste0("G:/scAPA/summary/Figure/fig4/tool.time.pptx"), width = 6, height = 4.5)

print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig4/tool.time.png"),width=6,height=4.5) 


##(1.2)memory usage-------------------------------------------------------------
p=ggplot(res, aes(x = `Cell Number`,y =`Max memory (MB)`/1024 , color =Method,fill= Method,group = Method)) +  #group=group：确保 geom_line() 按照 group 分组连接点
  #geom_line(linetype="dashed",size=1) +
  geom_line() +
  geom_point(size=1.8,shape=23,
             stroke = 0.5,            
             alpha=0.9) +  
  scale_color_manual(values=tool.col)+
  scale_fill_manual(values=tool.col)+
  labs(x = "Cell Number", y = "Max memory (GB)")+
  p_theme+
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1,size=7),
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    # 隐藏断点区域的 X 轴文本
    axis.text.x.top = element_blank(),     
    axis.ticks.x.top = element_blank(),     
    axis.line.x.top = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    legend.text = element_text(size = 7)
    
  )+
  #guides(fill="none")+
  scale_x_continuous(breaks = custom_breaks,  
                     labels = custom_breaks, 
                     limits = range(0,max(custom_breaks)))+  
  scale_y_continuous(breaks = seq(0,120,20), 
                     labels = seq(0,120,20), 
                     limits = range(0,120))+
  guides(color=guide_legend(title=''),fill='none')


topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig4/tool.mem.pptx"), width = 6, height = 4.5)


print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig4/tool.mem.png"),width=6,height=4.5) 


####2.Performance evaluation of pA identification-------------------------------------------------------


##（2.1）pA abundance score-------------------------------------------------------

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

#Summarized according to the sequencing protocol/species
group="protocol" #species/protocol



res <- sample.data %>%
  group_by(Species, Sample, Sample_ID) %>%  
  mutate(
    max.num = max(UTR3),
    min.num = min(UTR3),
    score = ifelse(
      max.num == min.num, 
      0.5,
      (UTR3 - min.num) / (max.num - min.num) 
    )
  ) %>%
  ungroup()

write.csv(res,file="G:/scAPA/summary/temp/fig4/pa.abundance.score.csv",row.names = FALSE)


if(group=="species"){
  mean.res<- res %>%
    group_by(Species, Tool) %>%  
    summarise(mean.score = mean(score, na.rm = TRUE)) %>% 
    ungroup() %>%
    arrange(Species,desc(mean.score))
}else if(group=="protocol"){
  mean.res<- res %>%
    group_by(Species, Protocol,Tool) %>%  
    summarise(mean.score = mean(score, na.rm = TRUE)) %>% 
    ungroup() %>%
    arrange(Species,Protocol,desc(mean.score))
}



#tool order-------------
all.mean.res= res %>% 
  group_by( Tool) %>% 
  summarise(all.mean.score = mean(score, na.rm = TRUE)) %>%  
  ungroup() %>%
  arrange(desc(all.mean.score))

#order=mean.res[mean.res$Species=="Human",]
order=all.mean.res

mean.res$Species<-factor(mean.res$Species,
                         levels = c("Human","Mouse","Arabidopsis")
)

mean.res$sub.group=paste0(mean.res$Species,"_",mean.res$Protocol)
mean.res$sub.group <- factor(mean.res$sub.group,
                             levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                             labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



mean.res$Tool=factor(mean.res$Tool,
                     levels = rev(order$Tool))
input.data=mean.res

###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab="pA Abundance Score"

p=ggplot(input.data, aes(x =sub.group, y =Tool, fill =mean.score)) + #fill =-1*log10(p.adj)
  geom_tile(
    color = "#FAFAFA",
    lwd = 0.8, 
    linetype = 1) + 
  #scale_fill_gradient2(low = '#3E4F94',
  #                     mid = 'white',
  #                     high = '#b5182b'
  #)+
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.score, 2))),  # 强制显示两位小数 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()       
        #axis.title.y = element_blank(),        
        #axis.text.y = element_blank(),       
        #axis.ticks.y = element_blank()       
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/",lab,"2.pptx"), width = 5.5, height = 5)


##(2.2)Precision and Overlap pA score----------------------------------------------------
#utr3.refpa <- readRDS("G:/scAPA/summary/temp/fig1/utr3.refpa(coord).rds")


pac.data=readRDS("G:/scAPA/summary/temp/fig1/pac.data_utr3.refpa(coord).rds")



pac.data$srr <- factor(pac.data$srr,
                       levels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                  "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                  "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                  "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                  "Arabidopsis1","Arabidopsis2","Arabidopsis3"))

pac.data$sub.species =paste0(pac.data$species,"_",pac.data$Protocol)

pac.data$sub.species=factor(pac.data$sub.species,
                            levels =c("Mouse_10x","Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq","Human_10x","Arabidopsis_10x" ),
                            labels = c("Mouse (10X)","Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)","Human (10X)","Arabidopsis (10X)"))


#pac.data$species <- factor(pac.data$species,
#                           levels=c("Human","Mouse","Arabidopsis"))

head(pac.data)

pac.data.need <- subset(pac.data,type=="Overlapping")
#pac.data.need$pro.value <- round(pac.data.need$Freq/pac.data.need$total*100,2)
#pac.data.need$rectheat <- 1

###The number of overlap pA segments standardized through min-max
pac.data.need <- pac.data.need %>%
  mutate(
    ref.pa.num = case_when(
      species == "Mouse"   ~ 132514,
      species == "Human"   ~ 117674,
      species == "Arabidopsis"  ~ 41267,
      TRUE                 ~ NA # 默认值
    )
  )
 
pac.data.need$recall=pac.data.need$Freq/pac.data.need$ref.pa.num
pac.data.need$precision=pac.data.need$Freq/pac.data.need$total


pac.data.need <- pac.data.need %>%
  group_by(species,srr) %>%  
  mutate(
    max.num = max(Freq),
    min.num = min(Freq),
    overlap.pa.score = ifelse(
      max.num == min.num,  # 如果组内数值全部相同，得分为0.5
      0.5,
      (Freq- min.num) / (max.num - min.num)  # min-max标准化
    )
  ) %>%
  ungroup()


index="Recall"

##tool order---------------
if(index=="Precision"){
  mean.res<- pac.data.need %>%
    group_by(species, Protocol, group) %>%  
    summarise(mean.precision = mean(precision, na.rm = TRUE)) %>%  
    ungroup() %>%
    arrange(species, Protocol,desc(mean.precision))
  
  order=pac.data.need %>% 
    group_by( group) %>%  
    summarise(mean.pro = mean(precision, na.rm = TRUE)) %>%  
    ungroup() %>%
    arrange(desc(mean.pro))
  
}else if(index=="Recall"){
  mean.res<- pac.data.need %>%
    group_by(species, Protocol, group) %>%  
    summarise(mean.score = mean(overlap.pa.score, na.rm = TRUE)) %>%  
    ungroup() %>%
    arrange(species, Protocol,desc(mean.score))
  
  order=pac.data.need %>% 
    group_by( group) %>%  
    summarise(mean.pro = mean(overlap.pa.score, na.rm = TRUE)) %>%  
    ungroup() %>%
    arrange(desc(mean.pro))
  
}


mean.res$sub.group=paste0(mean.res$species,"_",mean.res$Protocol)
mean.res$sub.group <- factor(mean.res$sub.group,
                             levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                             labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



mean.res$group =factor(mean.res$group,
                            levels =rev(order$group ))



###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab=index

p=ggplot(mean.res, aes(x =sub.group, y =group, fill =mean.score)) + #fill =-1*log10(p.adj)
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  #scale_fill_gradient2(low = '#3E4F94',
  #                     mid = 'white',
  #                     high = '#b5182b'
  #)+
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.score, 2))), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()       
       
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/pa.",lab,"2.pptx"), width = 5.5, height = 5)

write.csv(pac.data.need,file = "G:/scAPA/summary/temp/fig4/pac.data.need.csv",row.names = FALSE)


##(2.3)Accuracy of unique pA（0.5 PAS% +0.5 TP%）-----------------------------

##The score is calculated based on the proportion of polyA signals 
#and the proportion predicted as TP by DeepPASTA,with each accounting for 0.5 of the total score.
#The proportion of poly(A) signals------------------------
res=read.csv("G:/scAPA/summary/temp/fig1.2/uni.pa.pas.csv")

res <- res  %>% 
  dplyr::group_by(species,protocol,sample,sample_id,tool,PAS) %>% 
  dplyr::summarise(pas.num=sum(No),.groups = 'keep') 


res <- res %>% 
  dplyr::group_by(species,protocol,sample,sample_id,tool) %>% 
  dplyr::mutate(pas.fre=pas.num/sum(pas.num)) 

head(res)

# 定义必须包含的 PAS 类别
required_pas <- c("AATAAA", "1Variants", "Others")


res_complete <- as.data.frame(res) %>%
  tidyr::complete(
    tidyr::nesting(species, protocol, sample, sample_id, tool),
    PAS = required_pas,
    fill = list(pas.num = 0, pas.fre = 0)
  )

colnames(res_complete)[7] ="num"


data=res_complete %>%
  
  group_by(sample_id, tool) %>%
  
  mutate(
    pa.num = sum(num[1:min(3, n())]),
    
    
    pas.num = sum(num[PAS %in% c("AATAAA", "1Variants")])
  ) %>%
 
  mutate(pas.fre = sum(pas.fre[PAS %in% c("AATAAA", "1Variants")])) %>%
  select(-c(PAS,num)) %>%
  ungroup()

data=unique(data)


##TP proportion predicted by DeepPASTA------------------------

res=read.csv("G:/scAPA/summary/temp/fig1.2/deepPASTA.res_uni.pa.csv")

res=subset(res, !res$unique.pa.num==0 & res$prediction.label=="TRUE")

res$group=paste0(res$sample_id,"_",res$tool)


data$group=paste0(data$sample_id,"_",data$tool)


data <- data %>%
  left_join(res %>% select(group, fre), by = "group") %>%
  mutate(
    tp.fre = ifelse(is.na(fre), pas.fre, fre)  
  ) %>%
  select(-fre)

data$score=(data$pas.fre + data$tp.fre) /2

write.csv(data,file="G:/scAPA/summary/temp/fig4/unique.pa.acc.score.csv",row.names = FALSE)

###Calculate the average and sort------------------------------------
mean.res<- data %>%
  group_by(species, protocol,tool) %>%  
  summarise(mean.score = mean(score, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(species, protocol,desc(mean.score))


order=data %>% 
  group_by( tool) %>% 
  summarise(all.mean.score = mean(score, na.rm = TRUE)) %>%  
  ungroup() %>%
  arrange(desc(all.mean.score))


mean.res$sub.group=paste0(mean.res$species,"_",mean.res$protocol)
mean.res$sub.group <- factor(mean.res$sub.group,
                             levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                             labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



mean.res$tool=factor(mean.res$tool,
                     levels = rev(order$tool))



###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab="Unique pA\nAccuracy Score"

p=ggplot(mean.res, aes(x =sub.group, y =tool, fill =mean.score)) + #fill =-1*log10(p.adj)
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  #scale_fill_gradient2(low = '#3E4F94',
  #                     mid = 'white',
  #                     high = '#b5182b'
  #)+
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.score, 2))), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()       
        #axis.title.y = element_blank(),        
        #axis.text.y = element_blank(),       
        #axis.ticks.y = element_blank()      
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/unique.pa.accuracy.score",".pptx"), width = 5.5, height = 5)


###3.Performance of pA Quantification---------------------------------------------------------------
##（3.1）gene level:Correlation with gene count---------------------
res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/gcount.cor.csv")

mean.res=res %>%
  dplyr::group_by(species,protocol,method) %>%
  dplyr::summarise(
    mean.cor = mean(cor, na.rm = TRUE),
    #sd_value = sd(value, na.rm = TRUE),
    .groups = "keep"  
  )%>%
  arrange(species,protocol,method)

##tool order
order=mean.res %>% 
  group_by( method) %>%  
  summarise(all.mean.score = mean(mean.cor, na.rm = TRUE)) %>%  
  ungroup() %>%
  arrange(desc(all.mean.score))


mean.res$sub.group=paste0(mean.res$species,"_",mean.res$protocol)
mean.res$sub.group <- factor(mean.res$sub.group,
                             levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                             labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



mean.res$method=factor(mean.res$method,
                     levels = rev(order$method))



###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab="correlation coefficient"

p=ggplot(mean.res, aes(x =sub.group, y =method, fill =mean.cor)) + 
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.cor, 2))),   
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()        
        
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/gcount.cor2",".pptx"), width = 5.5, height = 5)
print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/gcount.cor",".png"),width=5.5,height=5)

##(3.2)Cell Clustering----------------------------------------------------------------
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.count.res.csv")

res=res[,c(1:4)]
res=subset(res,res$method %in% c("ARI","SC"))

res <- res %>%
  pivot_wider(
    names_from = method,  
    values_from = value   
  )


loc=match(res$sample,sample.data$Sample_ID)
res$protocol=sample.data$Protocol[loc]
res$species=sample.data$Species[loc]
write.csv(res,file="G:/scAPA/summary/temp/fig4/cluster.res.csv",row.names = FALSE)

index="SC"

mean.res=res %>%
  dplyr::group_by(species,protocol,tool) %>%
  dplyr::summarise(
    mean.ari = mean(ARI, na.rm = TRUE),
    mean.sc=mean(SC,na.rm=TRUE),
    #sd_value = sd(value, na.rm = TRUE),
    .groups = "keep"  
  )%>%
  arrange(species,protocol,tool)

##tool order
if(index=="ARI"){
  order=mean.res %>% 
    group_by( tool) %>%  
    summarise(all.mean.score = mean(mean.ari, na.rm = TRUE)) %>%  # 计算每组score的平均值
    ungroup() %>%
    arrange(desc(all.mean.score))
  
}else if(index=="SC"){
  
  order=mean.res %>% 
    group_by( tool) %>%  
    summarise(all.mean.score = mean(mean.sc, na.rm = TRUE)) %>%  # 计算每组score的平均值
    ungroup() %>%
    arrange(desc(all.mean.score)) 
  
}


mean.res$sub.group=paste0(mean.res$species,"_",mean.res$protocol)
mean.res$sub.group <- factor(mean.res$sub.group,
                             levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                             labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



mean.res$tool=factor(mean.res$tool,
                       levels = rev(order$tool))



###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab=index

p=ggplot(mean.res, aes(x =sub.group, y =tool, fill =mean.sc)) + 
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.sc, 2))), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()        
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/",index,".pptx"), width = 5.5, height = 5)
print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig4/score.heatmap/",index,".png"),width=5.5,height=5)

##(3.3)site level:Correlation of PAU----------------------------------
res=read.csv("E:/scAPA/summary/temp/fig2/bech_pacount/cor.res.csv")

#cor *0.3+ gene.num.score*0.7  
res <- res %>%
  group_by(sample) %>%  
  mutate(
    max.num = max(gene.num),
    min.num = min(gene.num),
    gene.num.score = ifelse(
      max.num == min.num,  
      0.5,
      (gene.num- min.num) / (max.num - min.num)  
    ),
    score = cor *0.3+ gene.num.score*0.7  
    
  ) %>%
  ungroup()%>%
  arrange(sample,desc(score))

write.csv(res,file="E:/scAPA/summary/temp/fig4/pau.cor.score.csv",row.names = FALSE)

###Calculate the average and tool order---------------------------------
mean.res<- res %>%
  group_by(sample,method) %>%  
  summarise(mean.score = mean(score, na.rm = TRUE)) %>%  
  ungroup() %>%
  arrange(sample,desc(mean.score))

order=res %>% 
  group_by( method) %>%  
  summarise(all.mean.score = mean(score, na.rm = TRUE)) %>%  
  ungroup() %>%
  arrange(desc(all.mean.score))

temp.res=order
temp.res$sample="Average"
temp.res=temp.res[,c(3,1,2)]
colnames(temp.res)[3]="mean.score"
mean.res=rbind(mean.res,temp.res)

mean.res$sample <- factor(mean.res$sample,
                             levels =c("mouse_sperm","mouse_tcell","mouse_tip","Average" ),
                             labels =c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Average") )



mean.res$method=factor(mean.res$method,
                     levels = rev(order$method))



###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab="Weighted correlation coefficient"

p=ggplot(mean.res, aes(x =sample, y =method, fill =mean.score)) + #fill =-1*log10(p.adj)
  geom_tile(
    color = "#FAFAFA",
    lwd = 0.8, 
    linetype = 1) + 
  #scale_fill_gradient2(low = '#3E4F94',
  #                     mid = 'white',
  #                     high = '#b5182b'
  #)+
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.score, 2))), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()       
  )

topptx(p, filename = paste0("E:/scAPA/summary/Figure/fig4/score.heatmap/pau.cor",".pptx"), width = 5.5, height = 5)
print(p)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig4/score.heatmap/pau.cor",".png"),width=5.5,height=5)

##4.Performance of the DEAPA gene identification----------------------------------------------------------

sample.data <- read.csv("E:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

res=read.csv("E:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.res.csv")

res$precision[is.na(res$precision)] <- 0

loc=match(res$sample,sample.data$Sample_ID)
res$protocol=sample.data$Protocol[loc]
res$species=sample.data$Species[loc]
write.csv(res,file="E:/scAPA/summary/temp/fig4/deapa.res.csv",row.names = FALSE)

##Calculate the average and tool order
index="Recall"

mean.res=res %>%
  dplyr::group_by(species,protocol,tool) %>%
  dplyr::summarise(
    mean.precision = mean(precision, na.rm = TRUE),
    mean.recall=mean(recall,na.rm=TRUE),
    #sd_value = sd(value, na.rm = TRUE),
    .groups = "keep"  
  )%>%
  arrange(species,protocol,tool)


if(index=="Precision"){
  order=res %>% 
    group_by(tool)%>%
    summarise(mean.precision=mean(precision, na.rm = TRUE),
              mean.recall=mean(recall, na.rm = TRUE))%>%
    arrange(desc(mean.precision))
} else if(index=="Recall"){
  order=res %>% 
    group_by(tool)%>%
    summarise(mean.precision=mean(precision, na.rm = TRUE),
              mean.recall=mean(recall, na.rm = TRUE))%>%
    arrange(desc(mean.recall))
}


mean.res$sub.group=paste0(mean.res$species,"_",mean.res$protocol)
mean.res$sub.group <- factor(mean.res$sub.group,
                             levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                             labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



mean.res$tool=factor(mean.res$tool,
                     levels = rev(order$tool))



###plot--------------------------------------------------------
col=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)

lab=index

p=ggplot(mean.res, aes(x =sub.group, y =tool, fill =mean.recall)) 
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  scale_fill_gradientn(
    colours = col(10))+  #col2(10)
  geom_text(aes(label = sprintf("%.2f", round(mean.recall, 2))),  
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill=lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank()        
  )

topptx(p, filename = paste0("E:/scAPA/summary/Figure/fig4/score.heatmap/de.",index,".pptx"), width = 4.5, height = 5)
print(p)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig4/score.heatmap/de.",index,".png"),width=4.5,height=5)


###5.Summary of various characteristics of the selected methods and their overall evaluation results----------------------
##based on the results of the analysis of the 10X Chromium datasets for human and mouse--------------------
path.lists=list()
path.lists[["pa.abundance.score"]]="E:/scAPA/summary/temp/fig4/pa.abundance.score.csv"
path.lists[["pa.accuracy"]]="E:/scAPA/summary/temp/fig4/pac.data.need.csv"
path.lists[["uniuqe.pa.accuracy"]]="E:/scAPA/summary/temp/fig4/unique.pa.acc.score.csv"

path.lists[["gcount.cor"]]="E:/scAPA/summary/temp/fig2/bech_pacount/gcount.cor.csv"
path.lists[["ari"]]="E:/scAPA/summary/temp/fig4/cluster.res.csv"
path.lists[["pau.cor"]]="E:/scAPA/summary/temp/fig4/pau.cor.score.csv"

path.lists[["de.accuracy"]]="E:/scAPA/summary/temp/fig4/deapa.res.csv"

indexs=names(path.lists)

res.lists=list()
mean.res.lists=list()

temp=data.frame(method=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe",
                         "MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"),
                Language=c("R","R","Python","Shell,R","Python,R",
                           "R","Shell,Python","R","Python,R,shell","Python,R"),
                Applicable.Data=c("All","All","All","Human/Mouse","All",
                                  "Human/Mouse","Human/Mouse","Human/Mouse","Human/Mouse","Human/Mouse")
)

for(i in 1:7){
  
  index=indexs[i]
  
  if(i==1){
    data=read.csv(path.lists[[index]])
    res=subset(data,data$Species %in% c("Mouse","Human")& data$Protocol=="10x")
    mean.res<- res %>%
      group_by(Tool) %>%  
      summarise(mean.score = mean(score, na.rm = TRUE)) %>%  # 计算每组score的平均值
      ungroup() %>%
      mutate(pa.abundance.rank = rank(-mean.score, ties.method = "min"))%>%
    arrange(pa.abundance.rank) 
    
    colnames(mean.res)[c(1,2)]=c("method",index)
    
    
  }else if(i==2){
    data=read.csv(path.lists[[index]])
    res=subset(data,data$species %in% c("Mouse","Human")& data$Protocol=="10x")
    
    mean.res<- res %>%
      group_by(group) %>%  
      summarise(real.pa.abundance.score=mean(overlap.pa.score,na.rm = TRUE),
        pa.precision = mean(precision, na.rm = TRUE )) %>%  
      ungroup() %>%
      #arrange(desc(precision))
    mutate(real.pa.abundance.rank = rank(- real.pa.abundance.score, ties.method = "min"),
           pa.precision.rank= rank(- pa.precision,ties.method = "min") )%>%
      arrange(real.pa.abundance.rank)
  
    colnames(mean.res)[1]="method"  
  }else if(i==3){
      
    data=read.csv(path.lists[[index]])
    res=subset(data,data$species %in% c("Mouse","Human")& data$protocol=="10x")
    
    mean.res<- res %>%
      group_by(tool) %>%  
      summarise(mean.score = mean(score, na.rm = TRUE)) %>% 
      ungroup() %>%
      #arrange(desc(mean.score))
      mutate(unique.pa.accuracy.rank = rank(- mean.score, ties.method = "min"))%>%
      arrange(unique.pa.accuracy.rank)
      
    colnames(mean.res)[c(1,2)]=c("method",index)
  }else if(i==4){
    data=read.csv(path.lists[[index]])
    res=subset(data,data$species %in% c("Mouse","Human")& data$protocol=="10x") 
    
    mean.res=res %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        gcount.cor = mean(cor, na.rm = TRUE),
        #sd_value = sd(value, na.rm = TRUE),
      )%>%
      ungroup() %>%
      mutate(gcount.cor.rank = rank(- gcount.cor, ties.method = "min"))%>%
      arrange(gcount.cor.rank)
    
  }else if(i==5){
      
    data=read.csv(path.lists[[index]])
    res=subset(data,data$species %in% c("Mouse","Human")& data$protocol=="10x")
    
    mean.res=res %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(
        ari = mean(ARI, na.rm = TRUE)
        #mean.sc=mean(SC,na.rm=TRUE),
        #sd_value = sd(value, na.rm = TRUE),
        
      )%>%
      ungroup() %>%
      mutate(ari.rank = rank(- ari, ties.method = "min"))%>%
      arrange(ari.rank)
    
    colnames(mean.res)[[1]]="method"
    
  }else if( i==6){
    data=read.csv(path.lists[[index]]) 
    
    mean.res<- data %>%
      group_by(method) %>%  
      summarise(pau.cor = mean(score, na.rm = TRUE)) %>%  
      ungroup() %>%
      mutate(pau.cor.rank = rank(- pau.cor, ties.method = "min"))%>%
      arrange(pau.cor.rank)
    
  }else if(i==7){
    data=read.csv(path.lists[[index]]) 
    res=subset(data,data$species %in% c("Mouse","Human")& data$protocol=="10x")
    
    mean.res=res %>%
      dplyr::group_by(tool) %>%
      dplyr::summarise(
        precision = mean(precision, na.rm = TRUE),
        recall=mean(recall,na.rm=TRUE),
        #sd_value = sd(value, na.rm = TRUE),
        #.groups = "keep"  
      )%>%
      ungroup() %>%
      mutate(precision.rank = rank(- precision, ties.method = "min"),
             recall.rank= rank(- recall,ties.method = "min") )%>%
      arrange(precision.rank)
    
    colnames(mean.res)[[1]]="method"
    }

    
    res.lists[[index]]=data   
    mean.res.lists[[index]]=mean.res
    
    temp=left_join(temp,mean.res,by="method") 
    
    message(paste0("Calculation ",index," ranking"))
}

saveRDS( res.lists,file="E:/scAPA/summary/temp/fig4/res.lists.rds")
saveRDS( mean.res.lists,file="E:/scAPA/summary/temp/fig4/mean.res.lists.rds")
write.csv(temp,file="E:/scAPA/summary/temp/fig4/summary.res.csv",row.names = FALSE)


##(5.2)plot------------------------------------------
library(funkyheatmap)


temp=read.csv("E:/scAPA/summary/temp/fig4/summary.res.csv")

order=temp %>%
  select(1:3,contains("rank")) %>%
  mutate(identification.rank = rowMeans(select(., 4:7), na.rm = TRUE),
         quantification.rank = rowMeans(select(., 8:10), na.rm = TRUE),
         de.rank = rowMeans(select(.,11:12), na.rm = TRUE)
        
         )%>%
rowwise() %>% 
  mutate(
    mean.rank = weighted.mean(
      c(identification.rank, quantification.rank, de.rank),
      c(0.4, 0.3, 0.3),
      na.rm = TRUE
    ),
   
    group = case_when(
      method %in% c("scAPAtrap", "Sierra", "SCAPE", "scAPA", "polyApipe") ~ "De-novo",
      method %in% c("MAAPER", "SCAPTURE", "Infernape", "scUTRquant", "scraps") ~ "Annotation-based",
      TRUE ~ NA_character_  # 默认值
    ))  %>%
  select( method,group,mean.rank)%>%
  group_by(group) %>%
  arrange(mean.rank, .by_group = TRUE) %>%
  ungroup()


temp=temp%>%
  arrange(match(method,order$method))

data=temp
colnames(data)[1]="id"
data$method=data$id
data=data %>%
  select(id,method,everything())

##Set group information for the ID
row_info<- data %>%
  mutate(
    id = id,  
    group = case_when(
      id %in% c("scAPAtrap", "Sierra", "SCAPE", "scAPA", "polyApipe") ~ "De-novo",
      id %in% c("MAAPER", "SCAPTURE", "Infernape", "scUTRquant", "scraps") ~ "Annotation-based",
      TRUE ~ NA_character_  
    )
  ) %>%
  select(id, group)


#column_info

need.col=data %>%
  select(-id, -contains("rank")) %>%
  colnames()

column_info=data.frame(id= need.col,
                       name=c("","Language","Species",
                              "pA abundance score","Abundance score of true pA","Precision","Uniuqe pA accuracy",
                              "Correlation with gene count","ARI","Correlation of PAU",
                              "Precision","Recall"),
                       group=c("name","name","name",
                               "identification","identification","identification","identification",
                               "count","count","count",
                               "DE", "DE"),
                       geom=c("text","text","text",
                              "funkyrect","funkyrect","funkyrect","funkyrect",
                              "funkyrect","funkyrect","funkyrect",
                              "funkyrect","funkyrect"),
                       palette=c(NA,"legend.col","legend.col",
                                 "Identification_Score","Identification_Score","Identification_Score","Identification_Score",
                                 "Quantification_Score","Quantification_Score","Quantification_Score",
                                 "DE_Score","DE_Score"),
                       
                       options = list(
                         list(legend = FALSE),
                         list(legend = TRUE),
                         list(legend = TRUE),
                         list(legend = FALSE), list(legend = FALSE), list(legend = FALSE), list(legend = FALSE),
                         list(legend = FALSE), list(legend = FALSE), list(legend = FALSE),
                         list(legend = FALSE), list(legend = FALSE)
                       )
)


col1=c("#fdeae0","#fad6c6","#f7bfa9","#f3a387","#f0896a")
col2=c("#f1fae0","#e4f1cc","#d2e3b7","#bad19c","#91af71")
col3=c("#c3edfc","#a4ddf8","#8bcbf1","#65b0e8","#498ac7")
col4=c("#edeff1","#d0d8dc","#afbec4","#91a4ae","#77919d")
col5=c("#d8f5e6","#bbe3d1","#9dc8b9","#77a69a","#558c83")
grey<- c("#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525")
#palettes <- list(
#  Identification_Score = grDevices::colorRampPalette(c(RColorBrewer::brewer.pal(9, "YlOrBr")[-7:-9]))(20),
#  Quantification_Score = grDevices::colorRampPalette(c(RColorBrewer::brewer.pal(9, "Greens")[-7:-9]))(20),
#  DE_Score = grDevices::colorRampPalette(c(RColorBrewer::brewer.pal(9, "Blues")[-7:-9]))(20)
#)

palettes <- list(
  Identification_Score = col1,
  Quantification_Score = col2,
  DE_Score = col3,
  legend.col=grey
)


##Legend Settings

legends <- list(
  list(
    geom = "funkyrect",
    palette = "legend.col",
    title = "Score"
    
  ),
  list(
    geom = "funkyrect",
    palette = "Quantification_Score",
    enabled=FALSE, 
    title = "",
    color="white"
    
  ),
  list(
    geom = "funkyrect",
    palette = "DE_Score",
    enabled=FALSE,
    title = "",
    color="white"
    
    
  ),
  list(
    geom = "funkyrect",
    palette = "Identification_Score",
    enabled=FALSE,
    title = "",
    color="white"
    
  )
  
)





column_groups=data.frame(group=c("name","name","name",
                                 "identification","identification","identification","identification",
                                 "count","count","count",
                                 "DE", "DE"),
                         Category=rep("\n", times =12 ),
                         Experiment=c("Method","Method","Method",
                                      "pA identification","pA identification","pA identification","pA identification",
                                      "pA quantification","pA quantification","pA quantification",
                                      "DEAPA gene\nidentification", "DEAPA gene\nidentification"),
                         
                         palette=c(NA,"legend.col","legend.col",
                                   "Identification_Score","Identification_Score","Identification_Score","Identification_Score",
                                   "Quantification_Score","Quantification_Score","Quantification_Score",
                                   "DE_Score","DE_Score")
)
p=funky_heatmap(
  data = data,
  column_info = column_info,
  column_groups = column_groups,
  #row_groups = row_groups,
  row_info = row_info,
  palettes = palettes
  ,
  legends=legends
  #position_args = position_arguments(col_annot_offset = 1.5)
)+ theme(legend.position = "none")+ 
  guides(fill = "none",color = "none")

topptx(p, filename = paste0("E:/scAPA/summary/Figure/fig4/score.heatmap/summary.res2",".pptx"), width = 12, height = 8)
print(p)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig4/score.heatmap/summary.res2",".png"),width=12,height=8)

