library(ggsci)
library(ggplot2)
library(eoffice)
#library(extrafont)
library(export)
library(rtracklayer)
library(GenomicFeatures)
library(movAPA)
library(wesanderson)
library(dplyr)
library(plyr)
library(patchwork)
library(aplot)
library(RColorBrewer)
library(pROC)
library(UpSetR)
library(psych)
library(ggpubr)
library(cowplot) 
library(Seurat)
library(cluster)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Athaliana.TAIR.TAIR9)




###plot parameter----------------------------------------------------------------------
p_theme<- theme_bw()+theme(
  text=element_text(family="Arial" ),
  axis.title.x =element_text(color="black",size=10,family="Arial" ) ,
  axis.title.y =element_text(color="black",size=10,family="Arial" ) ,
  axis.text.x =element_text(color="black",size=9,family="Arial" ) ,
  axis.text.y =element_text(color="black",size=9,family="Arial" ) ,
  legend.text =element_text(color="black",size=8,family="Arial"),
  legend.title=element_text(color="black",size=9,family= "Arial" ),
  legend.background = element_blank(),
  panel.border = element_blank(),panel.grid.major = element_blank(),
  panel.grid.minor =element_blank(),
  axis.line=element_line(colour = "black",size=0.4))


sample.col=c("Mouse_Sperm"="#775360","Mouse_T_Cell"="#e8b9c5","Mouse_TIP"="#a85559","Mouse_Stem_Cell"="#976E30","Mouse_GFP1"="#ecdfa6","Mouse_GFP2"="#ecdfa6",
             "Mouse_Intestine1"="#8080A3","Mouse_Intestine2"="#8080A3","Mouse_Intestine3"="#8080A3","Mouse_Intestine4"="#8080A3","Mouse_Intestine5"="#8080A3","Mouse_Intestine6"="#8080A3",
             "Mouse_Hypothalamus1"="#bfd3e1","Mouse_Hypothalamus2"="#bfd3e1","Mouse_Bone1"="#3d85b9","Mouse_Bone2"="#3d85b9",
             "Human_Pbmc1"="#133f7f","Human_Pbmc2"="#133f7f","Human_Pbmc3"="#133f7f","Human_Pbmc4"="#133f7f","Human_Covid_Pbmc1"="#CF6862","Human_Covid_Pbmc2"="#CF6862",
             "Arabidopsis1"="#577149","Arabidopsis2"="#577149","Arabidopsis3"="#577149")


#species.col=c("Mouse" ="#f9cb45","Human"='#b5182b',"Arabidopsis"="#387a88")

species.col=c("Mouse (10X)" ="#D94701","Mouse (CEL-seq)"="#FD8D3C","Mouse (Drop-seq)"="#FDBE85","Mouse (Microwell-seq)"="#ffd7c1","Human (10X)"='#163859',"Arabidopsis (10X)"="#064D4B") 

#####1.Analysis of Consensus pA-----------------------------------------------

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

GSElists=unique(sample.data$Sample)

data=data.frame()
new.paclists=list()
for(s in 1: length(GSElists)){
  
  sample=GSElists[s]
  species <- unique(sample.data[sample.data$Sample==sample,]$Species)
  #annoDB=anno.refpa[[species]]
  
  protocol <- unique(sample.data[sample.data$Sample==sample,]$Protocol)
  
  GSMlists=unique(sample.data[sample.data$Sample==sample,]$Sample_ID)
  
  ##不同样本3UTR PAC数据存放路径
  all.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  
  
  for(i in 1:length( GSMlists)) {
    #for(id in 15:15) {
    GSM <-  GSMlists[i]
    paclist=all.pac[[GSM]]

    # 遍历 paclist 中的每个元素
    for (j in 1:length(names(paclist))) {
     
      pac <- paclist[[j]]
      name=names(paclist)[j]
      
      new_counts <- matrix(1, 
                           nrow = nrow(pac@counts), 
                           ncol = 1,
                           dimnames = list(rownames(pac@counts),
                                           name))
      
                           
                           # 更新 counts
                           pac@counts <- new_counts
                           pac@colData=data.frame(row.names = name,group=name)
                           
                           # 更新回列表
                           paclist[[j]] <- pac
    }
    
    new.pac=mergePACds(paclist, d = 24, by = "coord")
    #head(new.pac@counts)
    #head(new.pac@anno)
    
    new.pac@anno$paid=paste0(new.pac@anno$chr,":",new.pac@anno$strand,":",new.pac@anno$coord)
    new.pac@counts=new.pac@counts[row.names(new.pac@anno),]
    row.names(new.pac@anno)=new.pac@anno$paid
    row.names(new.pac@counts)=row.names(new.pac@anno)
    
    new.paclists[[GSM]]=new.pac
    
    ###统计不同样本中pA位点被多少个工具所共同识别
    temp=data.frame(row.names =row.names(new.pac@counts),
                    species=species,
                    protocol=protocol,
                    sample=sample,
                    sample_id=GSM,
                    paid=row.names(new.pac@counts),
                    n_tools_detected =Matrix::rowSums(!new.pac@counts == 0)
                    )
    
    
    
    data=rbind(data,temp)
    
      
    message(paste0(sample,"(",GSM,")",": ",i,"/",length(GSMlists)))
  }
}

write.csv(data,file="G:/scAPA/summary/temp/fig1.2/com.pa.stat.csv",row.names = FALSE)
saveRDS(new.paclists,file="G:/scAPA/summary/temp/fig1.2/new.anno.paclists.rds")

head(new.paclists[[1]]@counts)
head(data)
table(data$n_tools_detected)

###plot---------------------------------
#data=read.csv("G:/scAPA/summary/temp/fig1.2/com.pa.stat.csv")
res<- data %>%
  group_by(species, protocol, sample, sample_id) %>%  
  dplyr::count(n_tools_detected, name = "count") %>%  
  mutate(fre = count / sum(count) ,
  percent = scales::percent(fre, accuracy = 0.01) 
 ) 

###Plot specified samples
res$sample_id <- factor(res$sample_id,
                                levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                         "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287",
                                         "GSM2333586","GSM2333587","GSM2906396","GSM2906399",
                                         "3k","4k","5k","6k","GSM4712885","GSM4712907",
                                         "GSM3490689","GSM3490690","GSM3490691"),
                                labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                           "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                           "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                           "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                           "Arabidopsis1","Arabidopsis2","Arabidopsis3"))

####(1)Grouped by species------------------------------------------------
cols=c("#60bfa4","#137f5e","#1e7a8d","#5b6e82","#c7d5b2","#e4d69f",
       "#ecbf51","#ed9d3a","#f06f25","#f73f43")
res$species=factor(res$species,
                   levels= c("Human","Mouse","Arabidopsis"))

res$sub.group=paste0(res$species,"_",res$protocol)
res$sub.group <- factor(res$sub.group,
                               levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                               labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )



p=ggplot(res, aes(x=n_tools_detected,y=round(fre,2),
                group=factor(n_tools_detected),
                color=factor(n_tools_detected))) +
  geom_boxplot(width=0.7,outlier.shape = 23,outlier.size = 1,outliers=F)+
  scale_color_manual(values=cols)+
  stat_boxplot(geom='errorbar',width=0.15)+
  labs(x="No. of tools",y=bquote("Percentage of pAs (%)"), #% of pAs
       title = "")+
  p_theme+
  scale_x_continuous(
    breaks = seq(1, 10, 1),
    labels = seq(1, 10, 1))+
  guides(fill="none",color="none")+
  theme(
    strip.text.x = element_text(size=9),
    element_rect(fill = "grey"))+
  facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)
  #facet_wrap(~ sub.group,scales="free_x",nrow=2,drop=TRUE)


topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/com.pa_species.pptx"), width =6, height = 4)
#保存为png
print(p)

graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/com.pa_species.png") ,width=6,height=4)

######（2）representative samples from different sequencing protocols-------------------------------
samples=c("Mouse_Sperm","Mouse_Bone1","Mouse_Hypothalamus1","Mouse_Intestine1","Human_Pbmc3","Arabidopsis3")
sub.res=subset(res,res$sample_id %in% samples)
table(sub.res$sample_id)




sub.res$sample_id <- factor(sub.res$sample_id,
                            levels=c("Mouse_Sperm","Human_Pbmc3","Arabidopsis3","Mouse_Intestine1","Mouse_Hypothalamus1","Mouse_Bone1"),
                            labels=c("Mouse_Sperm(10X)","Human_Pbmc3(10x)","Arabidopsis3(10x)",
                                     "Mouse_Intestine1(CEL-seq)","Mouse_Hypothalamus1(Drop-seq)","Mouse_Bone1(Microwell-seq)"))

###################################################################################################

#####（3）different sequencing protocols---------------------------------------------------
tech="Microwell-seq"
sub.res=subset(res,res$protocol %in% tech)
table(sub.res$sample_id)

############################################################################

head(sub.res)
cols=c("#60bfa4","#137f5e","#1e7a8d","#5b6e82","#c7d5b2","#e4d69f",
       "#ecbf51","#ed9d3a","#f06f25","#f73f43")

p=ggplot(data=sub.res, 
       aes(x=n_tools_detected, y=count,group=n_tools_detected,
           fill=factor(n_tools_detected))) +
  geom_bar(
    stat = "identity",
    width = 0.7,                  
    position = position_dodge(0.8) 
  )+
  
  geom_text(
    aes(label = count),  
    position = position_dodge(0.8), 
    vjust = -0.5,      
    size = 1.8         
  ) +
  p_theme+
  scale_x_continuous(
                     breaks = seq(1, 10, 1),
                     labels = seq(1, 10, 1))+ 
  labs(x="No. of tools",y=bquote("No. of pAs"),
       title = "")+
  #guides(fill='none')+  #隐藏图例
  guides(fill = guide_legend(title = 'No. of tools'))+
  scale_fill_manual(values=cols)+
  theme(
        strip.text.x = element_text(size=8),
        element_rect(fill = "grey")
        )+
  facet_wrap(~ sample_id,scales="free_x",nrow=2,drop=TRUE)

#topptx(p, filename = "G:/scAPA/summary/Figure/fig1.2/boxplot_com.pa/com.pa.sata.pptx", width = 8, height = 5)
topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_com.pa/com.pa_",tech,".pptx"), width = 8, height = 3)
#保存为png
print(p)
#graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_com.pa/com.pa.sata.png") ,width=8,height=5)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_com.pa/com.pa_",tech,".png") ,width=8,height=3)

####2.Agreement score----

#load data
data=read.csv("G:/scAPA/summary/temp/fig1.2/com.pa.stat.csv")

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

GSMlists=unique(sample.data$Sample_ID)

##（1）Calculate the agreement score for each sample---------------------------------
res=data.frame()
for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  GSM <-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
 
  temp=subset(data,data$sample_id==GSM)
  temp$pa.num=nrow(temp)
  temp$tools.num=length(sample.data[sample.data$Sample_ID==GSM,]$Tool)
 
  temp$tools.score=(temp$n_tools_detected-1)/(temp$tools.num-1)
  temp$score=sum(temp$tools.score)/temp$pa.num
  
  res=rbind(res,temp)
  
  message(paste0(sample,"(",GSM,")",": ",i,"/",length(GSMlists)))
}
  
write.csv(res,file="G:/scAPA/summary/temp/fig1.2/agreement.score.csv",row.names = FALSE)  

#res=read.csv("G:/scAPA/summary/temp/fig1.2/agreement.score.csv")
res=unique(res[,c(1:4,7,8,10)])

res$percent = scales::percent(res$score, accuracy = 0.01)


###（2）plot--------------------------------------------------------------
#seq.col=c("#775360",'#b5182b',"#a5cc26","#A46ACD","#E78D35","#002fa7")

res$sample_id <-factor(res$sample_id,
                       levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                "3k","4k","5k","6k","GSM4712885","GSM4712907","GSM3490689","GSM3490690","GSM3490691",
                                "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287",
                                "GSM2333586","GSM2333587",
                                "GSM2906396","GSM2906399"
                       ),
                       labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                  "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                  "Arabidopsis1","Arabidopsis2","Arabidopsis3",
                                  "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                  "Mouse_Hypothalamus1", "Mouse_Hypothalamus2",
                                  "Mouse_Bone1","Mouse_Bone2"))



res$protocol<-factor(res$protocol,
                     level=c("10x","Microwell-seq","Drop-seq","CEL-seq"),
                     labels = c("10X","Microwell-seq","Drop-seq","CEL-seq"))
seq.col=c("10X"="#a5cc26","CEL-seq"="#A46ACD","Drop-seq"="#E78D35","Microwell-seq"="#002fa7")


res$species<-factor(res$species,
                    level=c("Human","Mouse","Arabidopsis"))

test <- res %>%
  arrange(protocol,species,desc(score)) 


order=test$sample_id
res$sample_id= factor(res$sample_id,
                      levels =order)

#res$group=paste0(res$species,"_",res$protocol)

#res$group<-factor(res$group,
#                   level=c("Mouse_10x","Human_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"))


p=ggplot(data=res, 
       aes(x=sample_id, y=score,
           fill=protocol)) +
  geom_bar(
    stat = "identity"
    #width = 0.7,                
    #position = position_dodge(0.8)
  )+
  p_theme+
  labs(x="",y=bquote("Agreement Score "),
       title = "")+
  guides(fill = guide_legend(title = 'Protocol'))+
  #scale_fill_d3(palette = "category20")+
  scale_fill_manual(values=seq.col)+
  theme(
    axis.text.x = element_text(angle = 45,hjust=1,size = 7),
   
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/agreement.score2.pptx"), width = 7, height = 4)
#保存为png
print(p)

graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/agreement.score2.png") ,width=7,height=4)



####(2.1)Plot Legend------------------------------------------------
temp=res
temp$anno <- "Sample_ID"

p1.sample <-ggplot(  temp, aes(x=sample_id, y=anno, fill=sample_id)) + geom_tile() + 
  scale_y_discrete(position="left") +
  scale_fill_manual(breaks = temp$sample_id,
                    values=sample.col,
                    name="Sample_ID")+
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45,hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y=element_blank())+
  #axis.text.y =element_text(color="black",size=10,family="Helvetica" )) +
  xlab(NULL) + ylab(NULL) +guides(fill="none")




temp$anno2 <- "Species"


p1.special <- ggplot(temp, aes(x=sample_id, y=anno2, fill=species)) + geom_tile() + 
  scale_y_discrete(position="left") +
  scale_fill_manual(
    values=species.col,
    name="Species")+
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45,hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y =element_blank()) +
  xlab(NULL) + ylab(NULL) +guides(fill="none")

legend.data<-list(special=p1.special,one.sample=p1.sample)

sup.p<-(p + theme(axis.text.x = element_blank())) %>% 
  aplot::insert_bottom((legend.data$one.sample+theme(legend.position = "none")), height=.05)  %>% 
  aplot::insert_bottom((legend.data$special+theme(legend.position = "none")), height=.05)

topptx(sup.p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/sup.agreement.score.pptx"), width = 7, height = 4)
#保存为png
print(sup.p)

graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/sup.agreement.score.png") ,width=7,height=4)

############################################################################


####3.The number and proportion of unique pAs predicted by different methods----------------------------------
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
paclists=readRDS("G:/scAPA/summary/temp/fig1.2/new.anno.paclists.rds")

anno.refpa=readRDS("G:/scAPA/summary/temp/fig1/anno.utr3.refpa.rds")

GSMlists=names(paclists)


new.paclists2=list()
uni.paclists=list()
res=data.frame()
for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  GSM <-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)

  
  pac=paclists[[GSM]]
  
  temp <- data.frame(row.names = NULL,
    species=species,
    protocol=protocol,
    sample=sample,
    sample_id=GSM,
    tool = colnames(pac@counts),  
    pA.num_one.merge = colSums(pac@counts != 0)  
  )
  
  
  annodb=anno.refpa[[species]]
  
  # 创建新的 counts 矩阵
  counts <- matrix(1, 
                   nrow = nrow(annodb), 
                   ncol = 1,
                   dimnames = list(rownames(annodb),
                                       "refPA"))
  anno=createPACdataset(counts,annodb)
  
  ##合并来自不同工具的pA以及refPA
  temp.list=list()
  temp.list[["PA"]]=pac
  temp.list[["refPA"]]=anno
  
  new.pac=mergePACds(temp.list, d = 24, by = "coord")
  
  new.pac@anno$paid=paste0(new.pac@anno$chr,":",new.pac@anno$strand,":",new.pac@anno$coord)
  new.pac@counts=new.pac@counts[row.names(new.pac@anno),]
  row.names(new.pac@anno)=new.pac@anno$paid
  row.names(new.pac@counts)=row.names(new.pac@anno)
  new.paclists2[[GSM]]=new.pac
  
  data=new.pac@counts
  is.unique <- sapply(1:ncol(data), function(i) {
    (data[, i] > 0) & (Matrix::rowSums(data > 0) == 1)  
  })
  
 
  
  uni.count=NULL
  for(j in 1: (ncol(data)-1)){
    tool=colnames(data)[j]
    unique.site = row.names(data)[is.unique[, j]]
    sub.new.pac=subsetPACds(new.pac,PAs =unique.site )
    uni.paclists[[GSM]][[tool]]=sub.new.pac
    
    uni.count=c(uni.count,length(unique.site))
  }
  
  
  ###统计unique.PA数量
  temp$pA.num_twice.merge.by.refPA= colSums(data[,c(1:(ncol(data)-1))] != 0)
  temp$unique.count=uni.count
  
  res=rbind(res,temp)
  
  message(paste0(sample,"(",GSM,")",": ",i,"/",length(GSMlists)))
}

##计算unique pA的比例
res$unique.fre=res$unique.count/res$pA.num_twice.merge.by.refPA

write.csv(res,file="G:/scAPA/summary/temp/fig1.2/unique.pa.csv",row.names = FALSE)
 saveRDS(new.paclists2,file="G:/scAPA/summary/temp/fig1.2/new.anno.paclists2.rds")
saveRDS(uni.paclists,file="G:/scAPA/summary/temp/fig1.2/uni.paclists.rds")
  
###plot--------------------------------------------
 #res=read.csv("G:/scAPA/summary/temp/fig1.2/unique.pa.csv")
 
 
 res$overlap.count=res$pA.num_twice.merge.by.refPA-res$unique.count
 res$overlap.fre=1-res$unique.fre
 head(res)
 ##将数据转换为长格式，方便绘图
 input.data<- res %>%
   tidyr::pivot_longer(
     cols = c(unique.count, unique.fre, overlap.count, overlap.fre),
     names_to = c("type", ".value"),  
     names_sep = "\\." 
   ) %>%
   mutate(
     type = case_when(
       type == "unique" ~ "unique pA",
       type == "overlap" ~ "overlap pA"
     )
   )
 
 colnames(input.data)[7]="pA.num"
 
 ###Plot specified samples
 input.data$sample_id <- factor(input.data$sample_id,
                                levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                         "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287",
                                         "GSM2333586","GSM2333587","GSM2906396","GSM2906399",
                                         "3k","4k","5k","6k","GSM4712885","GSM4712907",
                                         "GSM3490689","GSM3490690","GSM3490691"),
                                labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                           "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                           "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                           "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                           "Arabidopsis1","Arabidopsis2","Arabidopsis3"))
#input.data$tool<- factor(input.data$tool,
#        levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
 
 
 
 ######（1）representative samples from different sequencing protocols-------------------------------
 samples=c("Mouse_Sperm","Mouse_Bone1","Mouse_Hypothalamus1","Mouse_Intestine1","Human_Pbmc3","Arabidopsis3")
 sub.input.data=subset(input.data,input.data$sample_id %in% samples)
 table(sub.input.data$sample_id)
 
 
 
 
 sub.input.data$sample_id <- factor(sub.input.data$sample_id,
                             levels=c("Human_Pbmc3","Mouse_Sperm","Arabidopsis3","Mouse_Intestine1","Mouse_Hypothalamus1","Mouse_Bone1"),
                             labels=c("Human_Pbmc3 (10X)","Mouse_Sperm (10X)","Arabidopsis3 (10X)",
                                      "Mouse_Intestine1 (CEL-seq)","Mouse_Hypothalamus1 (Drop-seq)","Mouse_Bone1 (Microwell-seq)"))
 
 
 ##tool order
 
 sub.input.data <- sub.input.data %>%
  
   arrange(sample_id, desc(pA.num)) %>%

   group_by(sample_id) %>%
   mutate(
    tool.order = factor(tool, levels = unique(tool))
   ) %>%
    ungroup()
   
 
 
 
 
 
 ###################################################################################################
 
 
 
 
 #####（2）different sequencing protocols---------------------------------------------------
 tech="10x"
 sub.input.data=subset(input.data,input.data$protocol %in% tech)
 table(sub.input.data$sample_id)
 
 ############################################################################
 
 sub.input.data$type <-factor(sub.input.data$type,
                        levels =rev(c("overlap pA","unique pA")),
                        labels = rev(c("Overlapping","Unique"))
                         )
 

  cols=c("Unique"="#CBC9E2","Overlapping"="#9392a9")
 
 text.data=subset(sub.input.data,sub.input.data$type=="Unique")
 text.data$y.label=text.data$pA.num-(text.data$count)*0.5
 
 p_theme= theme_bw()+theme(
   text=element_text(family="Arial" ),
   axis.title.x =element_text(color="black",size=10,family="Arial" ) ,
   axis.title.y =element_text(color="black",size=10,family="Arial" ) ,
   axis.text.x =element_text(color="black",size=8,family="Arial" ) ,
   axis.text.y =element_text(color="black",size=9,family="Arial" ) ,
   legend.text =element_text(color="black",size=8,family="Arial"),
   legend.title=element_text(color="black",size=9,family= "Arial" ),
   legend.background = element_blank(),
   panel.border = element_blank(),panel.grid.major = element_blank(),
   panel.grid.minor =element_blank(),
   axis.line=element_line(colour = "black",size=0.4))
 
 
 p=ggplot(data=sub.input.data,
        aes(x=tool.order,  y=count, fill=type)) +
   geom_bar(stat="identity",
            width = 0.7,                 
            position = "stack" 
            )+
   p_theme+
   labs(x=NULL,y="No. of pAs",title = "")+
   guides(fill=guide_legend(title=NULL))+
   scale_fill_manual(values=cols)+
  
   geom_text(data=text.data,aes(y=pA.num, label=paste0((round(fre , 2))* 100, "%")), 
             #position = position_stack(vjust = 0.5), 
             #vjust = -0.2,
             hjust = -0.1,
             color="#CBC9E2", 
             size=2)+
   theme(#legend.position = "top",
     axis.text.x = element_text(angle = 45,hjust=1),
     strip.text.x = element_text(size=8),
     element_rect(fill = "grey"))+
   #facet_wrap(~ sample_id,scales="free_x",nrow=2)
    coord_flip()+
   facet_wrap(~ sample_id,nrow=2)
 
 
 topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/uni.pa.sata2.pptx"), width = 8, height = 5)
 #topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/com.pa_",tech,".pptx"), width = 8, height = 12)
 #保存为png
 print(p)
 graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/uni.pa.sata2.png") ,width=8,height=5)
#graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/com.pa_",tech,".png") ,width=8,height=12)


###(3)The percentage of unique pAs containing poly(A) signal----------------------------------

setwd( "G:/scAPA/summary/temp/fig1.2/unique.pa_fa")
##注释unique pA

bsgenome.list <- list()
bsgenome.list[["Human"]] <- BSgenome.Hsapiens.UCSC.hg38
bsgenome.list[["Mouse"]] <- BSgenome.Mmusculus.UCSC.mm10
bsgenome.list[["Arabidopsis"]] <- BSgenome.Athaliana.TAIR.TAIR9


sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
paclists=readRDS("G:/scAPA/summary/temp/fig1.2/uni.paclists.rds")
annotation.lists<-readRDS("G:/scAPA/summary/temp/fig1/annotation.lists.rds")
GSMlists=names(paclists)


anno.paclists=list()
fafile.lists=list()
res=data.frame()
for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  GSM <-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  anno=annotation.lists[[species]]
  bsgenome=bsgenome.list[[species]]
  
  tool.lists=names(paclists[[GSM]])
  
  for(j in 1: length(tool.lists)){
    tool=tool.lists[j]
    pac=paclists[[GSM]][[tool]]
    
    if(nrow(pac@anno)==0){
      next
    }else{
      
      if(GSM=="GSM1524285" & tool=="scraps"){
        raw.pac=readRDS("G:/scAPA/summary/PACds/mouse_int_utr.paclists.rds")
        raw.pac=raw.pac[[GSM]][[tool]]
        raw.anno=raw.pac@anno
        row.names(raw.anno)=paste0(raw.anno$chr,":",raw.anno$strand,":",raw.anno$coord)
        raw.anno=raw.anno[row.names(pac@anno),]
        pac@anno=raw.anno
        
      }else{
        pac <- annotatePAC(pac, aGFF = anno)
      }
      
      
      if(species=="Human"){
        pac <- ext3UTRPACds(pac,ext3UTRlen = 2000)
      }
      if(species=="Mouse"){
        pac <- ext3UTRPACds(pac,ext3UTRlen = 2000)
      }
      if(species=="Arabidopsis"){
        pac <- ext3UTRPACds(pac,ext3UTRlen = 1000)
    
        pac@anno$chr=paste0("Chr",pac@anno$chr)
      }
      
      
      pac<- annotateByPAS(pac, bsgenome, grams='AATAAA', from=-50, to=1, label=NULL)
      pac <- annotateByPAS(pac, bsgenome, grams='V1', from=-50, to=1, label=NULL)
      pac@anno$pA.signal <- "Others"
      pac@anno$pA.signal[which(!is.na(pac@anno$V1_dist))] <- "1Variants"
      pac@anno$pA.signal[which(!is.na(pac@anno$AATAAA_dist))] <- "AATAAA"
      
      anno.paclists[[GSM]][[tool]]=pac
      
      ###统计polyA信号分布
      temp.pas<- as.data.frame( table(pac@anno$pA.signal ))
      colnames(temp.pas) <- c("PAS","No")
      
      temp.pas$pa.unm=nrow(pac@anno)
      temp.pas$tool=tool
      temp.pas$sample_id=GSM
      temp.pas$sample=sample
      temp.pas$protocol=protocol
      temp.pas$species=species
      
      res=rbind(res,temp.pas)
      
    }
    
    
    message(paste0(GSM,": ",j,"/",length(tool.lists)))
  }
}
    
write.csv(res,file="G:/scAPA/summary/temp/fig1.2/uni.pa.pas.csv",row.names = FALSE)
saveRDS(anno.paclists,file="G:/scAPA/summary/temp/fig1.2/anno.uni.paclists.rds")  

##############################################################################

##(3.1)Plotting poly(A) signal===========================================
#
#
#
#res<-read.csv("G:/scAPA/summary/temp/fig1.2/uni.pa.pas.csv")

#res$tool <- factor(res$tool,
#                    levels=rev(c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
#)

##All signals
#不同样本数据中，各工具识别的unique pA所含有的数量
res <- res  %>% 
  dplyr::group_by(species,protocol,sample,sample_id,tool,PAS) %>% 
  dplyr::summarise(pas.num=sum(No),.groups = 'keep') 

#不同样本数据中，各工具识别的 unique pA附近不同polya信号的比例
res <- res %>% 
  dplyr::group_by(species,protocol,sample,sample_id,tool) %>% 
  dplyr::mutate(pas.fre=pas.num/sum(pas.num)) 

head(res)


required_pas <- c("AATAAA", "1Variants", "Others")

# 确保每个分组包含所有 PAS 类别
res_complete <- as.data.frame(res) %>%
  tidyr::complete(
    tidyr::nesting(species, protocol, sample, sample_id, tool),
    PAS = required_pas,
    fill = list(pas.num = 0, pas.fre = 0)
  )


# 验证：检查每组的 pas.fre 总和
test=res_complete %>%
  dplyr::group_by(species, protocol, sample, sample_id, tool) %>%
  dplyr::summarise(
    total.freq = sum(pas.fre),
    .groups = 'keep'
  )





res_complete$PAS <- factor(res_complete$PAS,
                           levels=rev(c("AATAAA","1Variants","Others")) ,
                           labels=rev(c("AATAAA","1nt-variants","No signals")))


#=============getting average==================================================
###(1)Grouped by species------------------------------------------------------
res.mean  <- res_complete %>% 
  dplyr::group_by(species,tool,PAS) %>% 
  dplyr::summarise(pas.fre.mean=mean(pas.fre),.groups = "keep")




test=res.mean %>%
  dplyr::group_by(species, tool) %>%
  dplyr::summarise(
    total.mean.freq = sum(pas.fre.mean),
    .groups = 'keep'
  )


res.mean <- res.mean[order(res.mean$tool,
                           #poly.signal.all.pa.average$pa_level,
                           res.mean$PAS,
                           decreasing = T),]

res.mean$sub.group=paste0(res.mean$species,"_",res.mean$tool)






order<- res_complete[res_complete$PAS=="No signals" & !res_complete$species=="Arabidopsis",] %>% 
  dplyr::group_by(tool) %>% 
  dplyr::summarise(pas.fre.mean=mean(pas.fre),.groups = "keep")%>%
  ungroup()%>%
  arrange(desc(pas.fre.mean))


res.mean$tool= factor(res.mean$tool,
                      levels =order$tool,
                      labels = order$tool)


res.mean$species<- factor(res.mean$species,
                          levels= c("Human","Mouse","Arabidopsis"))


###(2)Grouping by species and sequencing technology---------------------------------------

res.mean  <- res %>% 
  dplyr::group_by(species,protocol,tool,PAS) %>% 
  dplyr::summarise(pas.fre.mean=mean(pas.fre),.groups = "keep") 


res.mean <- res.mean[order(res.mean$tool,
                           #poly.signal.all.pa.average$pa_level,
                           res.mean$PAS,
                           decreasing = T),]

res.mean$group=paste0(res.mean$species,"_",res.mean$protocol)
res.mean$sub.group=paste0(res.mean$group,"_",res.mean$tool)


res.mean$group<- factor(res.mean$group,
                          levels =c("Mouse_10x","Human_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"))

head(res.mean)




###plot------------------------


head(res.mean)
data <- ddply(res.mean, .(sub.group),
              transform, 
              label_ypos=cumsum(pas.fre.mean) - 0.5*pas.fre.mean)

sub.data=subset(data,data$PAS=="No signals")



p=ggplot(data,
         aes(x=tool, y=pas.fre.mean, fill=PAS)) +
  geom_bar(stat="identity")+p_theme+
  labs(x=NULL,y="Percentage of poly(A) signals (%)")+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=rev(c("#992813","#a24f47","#eccab7")))+
  #scale_fill_manual(values=rev(c("#F2AD00","#ECCBAE","#D6D6D6")))+
  geom_text(data=sub.data,
            aes(y=label_ypos, label=paste0(round(pas.fre.mean,2)*100,"%")), vjust=1, 
            color="white", size=2.5)+
  coord_flip()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0,vjust=0.6),
        strip.text.x = element_text(size=8),
        element_rect(fill = "grey"))+
  guides(fill = guide_legend(reverse=TRUE,title=""))+
  scale_y_continuous(labels = scales::percent_format())+
  facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)
  #facet_wrap(~ group,scales="free_x",nrow=2,drop=TRUE)

topptx(p, filename = "G:/scAPA/summary/Figure/fig1.2/uni.pa.pas(species).pptx", width = 7, height = 4)
#topptx(p, filename = "G:/scAPA/summary/Figure/fig1.2/uni.pa.pas(species.tech).pptx", width = 7, height = 7)

print(p)
graph2png(file="G:/scAPA/summary/Figure/fig1.2/uni.pa.pas(species).png",width=7,height=4)
#graph2png(file="G:/scAPA/summary/Figure/fig1.2/uni.pa.pas(species.tech).png",width=7,height=7)


####4.The percentage of positive pAs among the unique pAs were predicted using DeepPASTA.-------------------------------------


sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
paclists=readRDS("G:/scAPA/summary/temp/fig1.2/uni.paclists.rds")


GSMlists=names(paclists)[1:22]


new.paclists=list()
res=data.frame()
for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  GSM <-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  
  
  tool.lists=names(paclists[[GSM]])
  
  for(j in 1: length(tool.lists)){
    tool=tool.lists[j]
    pac=paclists[[GSM]][[tool]]
                      

    file.path=paste0("G:/scAPA/summary/temp/fig1.2/deepPASTA/result/",GSM,"_",tool,"_200nt.3UTR.prediction.res.txt")
   
    
    if(file.exists(file.path)==TRUE){
    
      data=read.delim(file.path,header = FALSE)
      colnames(data)=c("seq","predict.value")
      
      data <- data %>%
        mutate(
          paid = stringr::str_remove(seq, "^>") %>%  
            stringr::str_split("_") %>%             
           purrr::map_chr(~ paste(.x[1:3], collapse = "_") ) 
        )
      data$paid=gsub("_",":",data$paid)
          
        ##处理deepPASTA的预测结果，若大于0.5，则为真的pA,反之则为假
       loc=match(pac@anno$paid,data$paid) 
         #table(is.na(loc))
         
         if(any(is.na(loc))==TRUE){
           stop("deepPASTA预测的pA与unique pA不一致")
         }
         
         pac@anno$deepPASTA.pre.value=data$predict.value[loc]
         
         pac@anno=pac@anno %>%
           mutate(label=if_else(deepPASTA.pre.value>0.5, TRUE,FALSE))
         
         new.paclists[[GSM]][[tool]]=pac
         
         # 统计FALSE和TRUE的数量
         false_count <- sum(pac@anno$label == FALSE, na.rm = TRUE)
         true_count <- sum(pac@anno$label == TRUE, na.rm = TRUE)
         
        
         temp <- data.frame(row.names = NULL,
                            species=species,
                            protocol=protocol,
                            sample=sample,
                            sample_id=GSM,
                            tool =tool,
                            prediction.label=c("TRUE","FALSE"),
                            number=c( true_count,false_count),
                            unique.pa.num=nrow(pac@anno)) 
         
        
    } else{
      
      new.paclists[[GSM]][[tool]]=pac
      
      temp <- data.frame(row.names = NULL,
                         species=species,
                         protocol=protocol,
                         sample=sample,
                         sample_id=GSM,
                         tool =tool,
                         prediction.label=c("TRUE","FALSE"),
                         number=c(0,0),
                         unique.pa.num=nrow(pac@anno))
      
    }
    
    res=rbind(res,temp)
      
    message(paste0(GSM,": ",j,"/",length(tool.lists)))
  }
}                                 


  res$fre=res$number/res$unique.pa.num
  res$percent= scales::percent(res$fre, accuracy = 0.01) 


saveRDS(new.paclists,file="G:/scAPA/summary/temp/fig1.2/deepPASTA_uni.paclists.rds")
write.csv(res,file="G:/scAPA/summary/temp/fig1.2/deepPASTA.res_uni.pa.csv",row.names = FALSE)


####（4.1）plot-------------------------------------------------
res=read.csv("G:/scAPA/summary/temp/fig1.2/deepPASTA.res_uni.pa.csv")
res=subset(res, !res$unique.pa.num==0)

###Plot specified samples
res$sample_id <- factor(res$sample_id,
                               levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                        "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287",
                                        "GSM2333586","GSM2333587","GSM2906396","GSM2906399",
                                        "3k","4k","5k","6k","GSM4712885","GSM4712907"
                                        ),
                               labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                          "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                          "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                          "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2"
                                         ))

res$prediction.label<-factor(res$prediction.label,
                                 levels =c("FALSE","TRUE"),
                                 labels = c("Negative","Positive")
)

 
#####(1)Summarize by species (calculate the average）---------------------------------

res.mean  <- res %>% 
  dplyr::group_by(species,tool,prediction.label) %>% 
  dplyr::summarise(mean.fre=mean(fre),.groups = "keep") 


res.mean <- res.mean[order(res.mean$tool,
                           #poly.signal.all.pa.average$pa_level,
                           res.mean$prediction.label,
                           decreasing = T),]

res.mean$sub.group=paste0(res.mean$species,"_",res.mean$tool)


res.mean$species<- factor(res.mean$species,
                          levels= c("Human","Mouse","Arabidopsis"))


order=res[res$prediction.label=="Positive",] %>% 
  dplyr::group_by(tool) %>% 
  dplyr::summarise(mean.fre=mean(fre),.groups = "keep")%>%
  ungroup()%>%
  arrange(desc(mean.fre))


res.mean$tool= factor(res.mean$tool,
                      levels =rev(order$tool),
                      labels = rev(order$tool))




head(res.mean)
data <- ddply(res.mean, .(sub.group),
              transform, 
              #label_ypos=cumsum(mean.fre) - 0.5*mean.fre
              label_ypos=0.5*mean.fre)

##显示预测为真的比例
sub.data=subset(data,data$prediction.label=="Positive")



p=ggplot(data,
         aes(x=tool, y=mean.fre, fill=prediction.label)) +
  geom_bar(stat="identity")+p_theme+
  labs(x=NULL,y="Percentage of Positive pAs (%)")+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=c("Negative"="#BAE4B3","Positive"="#729a93"))+ 
  geom_text(data=sub.data,
            aes(y=label_ypos, label=paste0(round(mean.fre,2)*100,"%")), vjust=1, 
            color="white", size=2.5)+
  coord_flip()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0,vjust=0.6),
        #strip =FALSE, 
       
        strip.text.x = element_text(size=8),
        element_rect(fill = "grey"))+
  guides(fill = guide_legend(reverse=TRUE,title=""))+
  scale_y_continuous(labels = scales::percent_format())+
  facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)
#facet_wrap(~ group,scales="free_x",nrow=2,drop=TRUE)

topptx(p, filename = "G:/scAPA/summary/Figure/fig1.2/deepPASTA.res(species).pptx", width = 5, height = 4)

print(p)
graph2png(file="G:/scAPA/summary/Figure/fig1.2/deepPASTA.res(species).png",width=5,height=4)



######（2）representative samples from different sequencing protocols-------------------------------
samples=c("Mouse_Sperm","Mouse_Bone1","Mouse_Hypothalamus1","Mouse_Intestine1","Human_Pbmc3")
sub.res=subset(res,res$sample_id %in% samples)
table(sub.res$sample_id)




sub.res$sample_id <- factor(sub.res$sample_id,
                                   levels=c("Mouse_Sperm","Human_Pbmc3","Mouse_Intestine1","Mouse_Hypothalamus1","Mouse_Bone1"),
                                   labels=c("Mouse_Sperm(10X)","Human_Pbmc3(10x)",
                                            "Mouse_Intestine1(CEL-seq)","Mouse_Hypothalamus1(Drop-seq)","Mouse_Bone1(Microwell-seq)"))


###################################################################################################

#####（3）different sequencing protocols---------------------------------------------------
tech="Microwell-seq"
sub.res=subset(res,res$protocol %in% tech)
table(sub.res$sample_id)

###############################################################################################


cols=c("Negative"="#729a93","Positive"="#BAE4B3")


text.data=subset(sub.res,sub.res$prediction.label=="Positive")
#text.data$y.label=text.data$pA.num-(text.data$count)*0.5


test= res %>%
  arrange(sample_id, desc(unique.pa.num))

order=unique(test[test$sample_id == "Mouse_Sperm(10X)" | test$sample_id == "Mouse_Sperm" ,]$tool)

sub.res$tool=factor(sub.res$tool,
                    level=c(order,"scraps"),
                    labels = c(order,"scraps"))
                    


p_theme= theme_bw()+theme(
  text=element_text(family="Arial" ),
  axis.title.x =element_text(color="black",size=10,family="Arial" ) ,
  axis.title.y =element_text(color="black",size=10,family="Arial" ) ,
  axis.text.x =element_text(color="black",size=8,family="Arial" ) ,
  axis.text.y =element_text(color="black",size=9,family="Arial" ) ,
  legend.text =element_text(color="black",size=8,family="Arial"),
  legend.title=element_text(color="black",size=9,family= "Arial" ),
  legend.background = element_blank(),
  panel.border = element_blank(),panel.grid.major = element_blank(),
  panel.grid.minor =element_blank(),
  axis.line=element_line(colour = "black",size=0.4))


p=ggplot(data=sub.res,
         aes(x=tool, y=number, fill=prediction.label)) +
  geom_bar(stat="identity",
           width = 0.7,                 
           position = "stack" 
  )+
  p_theme+
  labs(x=NULL,y="No. of Unique pAs",title = "")+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=cols)+
  #设置文本标签位置
  geom_text(data=text.data,aes(y=unique.pa.num, label=paste0((round(fre , 2))* 100, "%")), 
            vjust = -0.2,
            color="#BAE4B3", 
            size=2.5)+
  theme(#legend.position = "top",
    axis.text.x = element_text(angle = 45,hjust=1),
    strip.text.x = element_text(size=8),
    element_rect(fill = "grey"))+
  facet_wrap(~ sample_id,scales="free_x",nrow=1,drop=TRUE)

#topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/deepPASTA.res.pptx"), width = 8, height = 5.5)
topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/deepPASTA.res_",tech,".pptx"), width = 8, height = 3)
#保存为png
print(p)
#graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/deepPASTA.res.png") ,width=8,height=5.5)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig1.2/boxplot_uni.pa/deepPASTA.res_",tech,".png") ,width=8,height=3)


