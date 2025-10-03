

rm(list=ls())
options(stringsAsFactors = FALSE)
#library(liftOver)
library(rtracklayer)
library(Biostrings)
library(plyr)
library(ggplot2)
library(eoffice)
#library(extrafont)
library(export)
library(rtracklayer)
library(GenomicFeatures)
library(movAPA)
library(wesanderson)
library(dplyr)
library(patchwork)
library(aplot)
library(RColorBrewer)
library(pROC)
library(UpSetR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Athaliana.TAIR.TAIR9)


source("G:/scAPA/code/function.R")
output <- "G:/scAPA/summary/Figure/fig1"
if(!dir.exists(output)){
  dir.create(output)
}

setwd(output)

###plot parameter-----------------------------------------------------------------
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

"#895c46"

sample.col=c("Mouse_Sperm"="#775360","Mouse_T_Cell"="#e8b9c5","Mouse_TIP"="#a85559","Mouse_Stem_Cell"="#976E30","Mouse_GFP1"="#ecdfa6","Mouse_GFP2"="#ecdfa6",
             "Mouse_Intestine1"="#8080A3","Mouse_Intestine2"="#8080A3","Mouse_Intestine3"="#8080A3","Mouse_Intestine4"="#8080A3","Mouse_Intestine5"="#8080A3","Mouse_Intestine6"="#8080A3",
             "Mouse_Hypothalamus1"="#bfd3e1","Mouse_Hypothalamus2"="#bfd3e1","Mouse_Bone1"="#3d85b9","Mouse_Bone2"="#3d85b9",
             "Human_Pbmc1"="#133f7f","Human_Pbmc2"="#133f7f","Human_Pbmc3"="#133f7f","Human_Pbmc4"="#133f7f","Human_Covid_Pbmc1"="#CF6862","Human_Covid_Pbmc2"="#CF6862",
             "Arabidopsis1"="#577149","Arabidopsis2"="#577149","Arabidopsis3"="#577149")


 
#species.col=c("Mouse" ="#f9cb45","Human"='#b5182b',"Arabidopsis"="#387a88")

species.col=c("Mouse (10X)" ="#D94701","Mouse (CEL-seq)"="#FD8D3C","Mouse (Drop-seq)"="#FDBE85","Mouse (Microwell-seq)"="#ffd7c1","Human (10X)"='#163859',"Arabidopsis (10X)"="#064D4B") 
species.col2=c("Mouse" ="#D99267","Human"='#163859',"Arabidopsis"="#064D4B")
species.col3=c("Mouse (10X)" ="#80546b","Mouse (CEL-seq)"="#ab7090","Mouse (Drop-seq)"="#d4a2bf","Mouse (Microwell-seq)"="#e2c7d8","Human (10X)"='#163859',"Arabidopsis (10X)"="#064D4B") 


tool.col=c("scAPAtrap"="#064D4B","Sierra"="#fbd26a","SCAPE"="#ABC8E5","scAPA"="#6666CC","polyApipe"="#e8b9c5","MAAPER"="#64d9b8",
           "SCAPTURE"="#b5182b","Infernape"="#775360","scUTRquant"="#7fb80e","scraps"="#FF6B00")

tool.col2=c("#064D4B","#fbd26a","#b5182b","#775360","#ABC8E5","#6666CC","#64d9b8","#e8b9c5","#7fb80e","#FF6B00")

tool.col3=c("#064D4B","#fbd26a","#002fa7","#6666CC","#b56387","#64d9b8","#b5182b","#775360","#7fb80e","#FF6B00")


heatmap.col=c("#313772","#2c4ca0","#326db6","#478ecc","#75b5dc",
              "#fee3ce","#eabaa1","#dc917b","#d16d5b","#c44438","#b7282e")

feature.col=c("#992813","#627b7f","#c38458","#9a8870","#caad9d")
base.col=c('#b5182b',"#387a88","#f9cb45","#133f7f")


##Loading Data----------------------------------------------------------------
#================================================================================
#
#       [1] Loading Data
#
#===============================================================================
#Sample_ID=c("GSM2803334+GSM280333","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848")
samples=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
          "arab")


all_GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
               "GSM2333586","GSM2333587", "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490690","GSM3490691")

#load data
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

table(sample.data$Tool)

#Infernape     MAAPER  polyApipe      scAPA  scAPAtrap      SCAPE   SCAPTURE     scraps scUTRquant     Sierra 
#22         25         19         24         25         24         22         25         12         25 


#
#
#====1.1:Order data from high to low based on selected tools

sample.data$Sample_ID <- factor(sample.data$Sample_ID,
                               levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                        "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","GSM2333586","GSM2333587","GSM2906396","GSM2906399",
                                        "3k","4k","5k","6k","GSM4712885","GSM4712907",
                                        "GSM3490689","GSM3490690","GSM3490691"),
                               labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                          "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                          "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                          "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                          "Arabidopsis1","Arabidopsis2","Arabidopsis3"))
#sample.data$Tool <- factor(sample.data$Tool,
#                           levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps")
#                            )

sample.data.sub <- subset(sample.data,Tool=="scAPAtrap")

sample.data.sub <- sample.data.sub[order(sample.data.sub$Species,
                                         -sample.data.sub$Total),]
sample.data.sub$Sample_ID <- factor(sample.data.sub$Sample_ID,
                                    levels= c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                              "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                              "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                              "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                              "Arabidopsis1","Arabidopsis2","Arabidopsis3"))

##Plotting Legend--------------------------------------------------------
#================================================================================
#
#       [2] Plotting Legend
#
#===============================================================================

#for keeping consistent color
#sample.info.sub.pac5 <- sample.info.sub[order(sample.info.sub$Species,
#                                              -sample.info.sub$PAC_5),]

sample.data.sub2 <- subset(sample.data,Tool=="scAPAtrap")
sample.data.sub2$Sample_ID <- factor(sample.data.sub2$Sample_ID,
                                     levels= c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                               "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                               "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                               "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                               "Arabidopsis1","Arabidopsis2","Arabidopsis3"))
legend.data <- sample.data[!duplicated(sample.data$Sample_ID),]
#legend.data  <- legend.data[order(legend.data$Protocol),]
legend.data$Sample_ID <- factor(legend.data$Sample_ID,
                                levels= c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                          "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                          "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                          "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                          "Arabidopsis1","Arabidopsis2","Arabidopsis3"))

legend.data$sub.sample=paste0(legend.data$Species,"_",legend.data$Protocol)
legend.data$sub.sample= factor(legend.data$sub.sample,
                               levels =c("Mouse_10x","Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq","Human_10x","Arabidopsis_10x" ),
                               labels = c("Mouse (10X)","Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)","Human (10X)","Arabidopsis (10X)"))

legend.data$Species <- factor(legend.data$Species,
                                levels= c("Mouse","Human","Arabidopsis"))

legend.data$anno <- "Sample_ID"
one.sample <-ggplot(legend.data, aes(Sample_ID, y=anno, fill=Sample_ID)) + geom_tile() + 
  scale_y_discrete(position="left") +
  #scale_fill_manual(breaks = (sample.info.sub.pac5$Protocol),
  #                  values=c("#293352","#446455","#009999", "#972D15",
  #                           RColorBrewer::brewer.pal(n = 12, name = "Set3"),
  #                           "#F7F7F7","#67A9CF"),
   #                 name="Protocal")+
  
  scale_fill_manual(breaks = levels(sample.data.sub2$Sample_ID),
                    values=sample.col,
                   name="Sample_ID")+
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45,hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y =element_text(color="black",size=9,family="Arial" )) +
  xlab(NULL) + ylab(NULL) #+guides(fill=guide_legend(nrow=2,byrow=TRUE))

legend.data$anno2 <- "Species"
one.special <-ggplot(legend.data, aes(Sample_ID, y=anno2, fill=Species)) + geom_tile() + 
  scale_y_discrete(position="left") +
  #scale_fill_viridis_d(name = "Species", 
  #                    labels = c("Mouse","Human",  "Arabidopsis"), 
  #                     option = "viridis")+
  scale_fill_manual(
                    values=species.col2,
                    name="Species")+
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45,hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y =element_text(color="black",size=9,family="Arial" )) +
  xlab(NULL) + ylab(NULL) 

legend.data$anno3 <- "sub.sample"
sub.sample <-ggplot(legend.data, aes(Sample_ID, y=anno2, fill=sub.sample)) + geom_tile() + 
  scale_y_discrete(position="left") +
  #scale_fill_viridis_d(name = "Species", 
  #                    labels = c("Mouse","Human",  "Arabidopsis"), 
  #                     option = "viridis")+
  scale_fill_manual(
    values=species.col,
    name="Species")+
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45,hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y =element_text(color="black",size=9,family="Arial" )) +
  xlab(NULL) + ylab(NULL) 

legend.data<-list(special=one.special,sub.sample=sub.sample,one.sample=one.sample)

#===plot
leg.sample <- ggpubr::get_legend(legend.data$one.sample+
                                     theme(legend.justification = "left")+
                                     guides(fill=guide_legend(nrow=5,byrow=TRUE)))
ggpubr::as_ggplot(leg.sample)

leg.sub.sample <- ggpubr::get_legend(legend.data$sub.sample+
                                    theme(legend.justification = "left")+
                                    #theme(legend.position = c(1, 1))+
                                    guides(fill=guide_legend(nrow=1,byrow=TRUE)))
ggpubr::as_ggplot(leg.sub.sample)

leg.species <- ggpubr::get_legend(legend.data$special+
                                    theme(legend.justification = "left")+
                                    #theme(legend.position = c(1, 1))+
                                    guides(fill=guide_legend(nrow=1,byrow=TRUE)))
ggpubr::as_ggplot(leg.species)

aplot::plot_list(ggpubr::as_ggplot(leg.sub.sample) ,ggpubr::as_ggplot(leg.sample)  ,
                 nrow=2, tag_levels= NULL)
graph2png(file="legend_bar_annotation.png",
          width=9,height=3)

p=aplot::plot_list(ggpubr::as_ggplot(leg.sub.sample) ,ggpubr::as_ggplot(leg.sample)  ,
                 nrow=2, tag_levels= NULL)
topptx(p, filename ="legend_bar_annotation.pptx", width = 9, height = 3,units = "in")

##plot fig1a:number of pac-------------------------------------------------------
#================================================================================
#
#       [3] Plotting Figure1.A: number of pac
#
#===============================================================================

###根据识别的pA总量，指定工具顺序
data <- sample.data %>%
  group_by(Tool) %>%
  mutate(total_pa = sum(Total),
    total_utr3_pa = sum(UTR3)) %>%
  ungroup() %>%
  arrange(desc(total_pa))



sample.data$Tool=factor(sample.data$Tool,
                        levels = rev(unique(data$Tool)),
                        labels =rev(unique(data$Tool))
                        )



p1<- ggplot(data=sample.data, 
                aes(x=Sample_ID, y=Total/1000, group=Sample_ID,fill=Sample_ID)) +
  geom_bar(stat="identity")+
  p_theme+
  labs(x=NULL,y=bquote("# pA (x 10"^3~")"))+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  guides(fill='none')+
  scale_fill_manual(values=sample.col)+
  facet_grid(vars(Tool),scales="free_x",drop=TRUE,space="free")+
  theme(strip.text.y = element_text(size = 7))

topptx(p1, filename ="fig1.pac_number.pptx", width = 9, height = 7.5,units = "in")
##保存为png
print(p1)
graph2png(file="fig1.pac_number.png",width=9,height=7.5)


sup.p1<-(p1 + theme(axis.text.x = element_blank())) %>% 
  aplot::insert_bottom((legend.data$one.sample+theme(legend.position = "none")), height=.05)  %>% 
  aplot::insert_bottom((legend.data$sub.sample+theme(legend.position = "none")), height=.05)

topptx(sup.p1, filename ="fig1.sup.pac_number.pptx", width = 9, height = 8,units = "in")

print(sup.p1)
graph2png(file="fig1.sup.pac_number.png",width=9,height=8)



##plot fig1b:number of 3UTR.PAC---------------------------------------------
#================================================================================
#
#       [4] Plotting Figure1.B: number of 3UTR.PAC
#
#===============================================================================

data <- sample.data %>%
  group_by(Tool) %>%
  mutate(total_pa = sum(Total),
         total_utr3_pa = sum(UTR3)) %>%
  ungroup() %>%
  arrange(desc(total_utr3_pa))



sample.data$Tool=factor(sample.data$Tool,
                        levels = rev(unique(data$Tool)),
                        labels =rev(unique(data$Tool))
                        )



p2 <- ggplot(data=sample.data, 
                 aes(x=Sample_ID, y=UTR3/1000, group=Sample_ID,fill=Sample_ID)) +
  geom_bar(stat="identity")+
  p_theme+
  labs(x=NULL,y=bquote("# 3UTR pA (x 10"^3~")"))+
  theme(axis.text.x = element_text(angle = 45,hjust=1),
        strip.text.y = element_text(size = 7.5))+
  #guides(fill=guide_legend(title=NULL))+
  guides(fill='none')+
  scale_fill_manual(values=sample.col)+
  facet_grid(vars(Tool),scales="free_x",drop=TRUE,space="free")

topptx(p2, filename ="fig1.utr3.pac_number.pptx", width = 9, height = 7.5,units = "in")
##保存为png
print(p2)
graph2png(file="fig1.utr3.pac_number.png",width=9,height=7.5)


sup.p2 <-(p2 + theme(axis.text.x = element_blank())) %>% 
  insert_bottom((legend.data$one.sample+theme(legend.position = "none")), height=.05)  %>% 
  insert_bottom((legend.data$sub.sample+theme(legend.position = "none")), height=.05)

topptx(sup.p2, filename ="fig1.sup.utr3.pac_number.pptx", width = 9, height = 8,units = "in")

print(sup.p2)
graph2png(file="fig1.sup.utr3.pac_number.png",width=9,height=8)


####plot 1c: pA Abundance Score----------------------------------------------
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)


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

p=ggplot(data=mean.res, 
       aes(x=mean.score, y= Tool,group=Tool,fill=Tool)) +
  geom_bar(stat="identity",alpha=0.7)+
  p_theme+
  
  geom_text(data=mean.res,aes(y= Tool, label=round(mean.score , 2)), 
            #position = position_stack(vjust = 0.5),  # 堆叠居中
            hjust =  -0.2,
            color="black", 
            size=2.5)+
  labs(y=NULL,x="pA Abundance Score")+
  theme(#axis.text.x = element_text(angle = 45,hjust=1),
        strip.text.y = element_text(size = 7.5))+
  guides(fill='none')+
  scale_fill_manual(values=tool.col)+
  facet_wrap(~ sub.group,nrow=2,drop=TRUE)           


topptx(p, filename =paste0( "fig1.pa.num.score(",group,")",".pptx"), width = 7, height = 8)

#保存为png
print(p)
graph2png(file=paste0( "fig1.pa.num.score(",group,")",".png") ,width=7,height=8)


##plot2.Density plot between predicted PAC and known PACs-------------------------------------------
#===Loading All predicted PAS and PAC list ====================================
#

tool.list_all <- readRDS("G:/scAPA/summary/temp/fig1/gr_paclists(coord).rds")
tool.list <- readRDS("G:/scAPA/summary/temp/fig1/utr3_gr_paclists(coord).rds")

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

##this is collected reference PAC table
refpa <- readRDS("E:/scAPA/summary/temp/fig1/refpa(coord).rds")
utr3.refpa <- readRDS("E:/scAPA/summary/temp/fig1/utr3.refpa(coord).rds")


#The nearest distance between the 3'UTR PA site identification by each tool and the reference PA
samples=unique(sample.data$Sample_ID)

pac.distance <- data.frame()
new.tool.list=list()
for( id in 1:length(samples)){
  sample_id <- samples[id]
  species <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Protocol)
  tool <- tool.list[[sample_id]]
  for(tools in names(tool)){
    
    ###all pac
    pac_all <- tool.list_all[[sample_id]][[tools]]
    pac.result_all <- getDistance.region(query.pa=pac_all,
                                     known.pa=refpa[[species]],
                                     srrID=sample_id,
                                     group=tools,
                                     species=species,
                                     protocol=protocol)
     new.tool.list[[sample_id]][[tools]]$pac <- pac.result_all$query.data 
    
    # delete pas table if you don't need it anymore
    #tool.list[[srrID]][[tools]]$pas <- NULL
    
    ###utr3 pac
    pac <- tool.list[[sample_id]][[tools]]
    pac.result <- getDistance.region(query.pa=pac,
                                     known.pa=utr3.refpa[[species]],
                                     srrID=sample_id,
                                     group=tools,
                                     species=species,
                                     protocol=protocol)
    new.tool.list[[sample_id]][[tools]]$utr3.pac <- pac.result$query.data   
    
    pac.distance <-rbind(pac.distance,pac.result$distance.data)
    
    print(paste0(sample_id,": ",id,"/",length(samples)," (run ",tools," pac data)"))
  }
}

saveRDS(new.tool.list,file="G:/scAPA/summary/temp/fig1/gr_paclists_addDistance(coord).rds")

saveRDS(pac.distance,file="G:/scAPA/summary/temp/fig1/get_pac.distance(coord).rds")

###plot density plot-------------------------------------------------------------------------
#pac.distance=readRDS("E:/scAPA/summary/temp/fig1/get_pac.distance(coord).rds")
#dim(pac.distance)
#[1] 2565659       6



#pac.distance$group <- factor(pac.distance$group,
#                             levels=c("MAAPER","SCAPTURE","Infernape","scUTRquant","scraps","scAPAtrap","Sierra","SCAPE","scAPA","polyApipe"))

##Grouping by species----------

input.data <- pac.distance
input.data  <- input.data[,c(1,2,4,5,6)]
#input.level="PAC3"
#cutoff.count <- 3
#input.data <- subset(input.data,count>=cutoff.count)
input.data$count <- 1


#tool order
order=input.data %>%
  group_by(species,group,distance)%>%
 summarise(count=sum(count))%>%
  ungroup()%>%
  group_by(group) %>%
 
  summarise(
    total.count=sum(count),
    count_in_range = sum(count[distance >= -50 & distance <= 50], na.rm = TRUE),
    
    fre = count_in_range / total.count,
    .groups = "drop"
  ) %>%
  arrange(desc(fre)) 

              


input.data <- input.data %>%
  #dplyr::group_by(species,protocol,group,distance) %>%
  dplyr::group_by(species,group,distance) %>%
  dplyr::summarise(count=sum(count))

#input.data$sub.group=paste0(input.data$species,"_",input.data$protocol)


input.data <- input.data %>%
  #dplyr::group_by(species,protocol,group) %>%
  dplyr::group_by(species,group) %>%
  dplyr::mutate(total.count=sum(count))

input.data$fren <-   input.data$count/input.data$total.count

input.data <-   input.data[order(  input.data$fren,decreasing = T),]



input.data$species <- factor(input.data$species,
                          levels =c("Human","Mouse","Arabidopsis"))
input.data$group <- factor(input.data$group,
                           levels = order$group)


#near=subset(input.data,abs(input.data$distance)<=100)
#near_per=round(sum(near$count)/unique(near$total.count),4)
#scales::percent(near_per)/sprintf("%1.2f%%", 100*near_per)

p<- ggplot(input.data, aes(x=distance,weight=fren, 
                           fill=group,
                           color=group)) +
  geom_density(size=0.8,alpha=0.2)+scale_x_continuous(limits = c(-200, 200))+
  p_theme+
  scale_fill_manual(values=tool.col)+
  scale_color_manual(values=tool.col)+
  labs(x="Distance (nt)",y="Density")+
  theme(#legend.position  =c(0.2,0.7),
        axis.text.x = element_text(angle = 0,vjust=0.6),
        #strip =FALSE, 
        strip.text.x = element_text(size=9),
        element_rect(fill = "grey"))+
  guides(color=guide_legend(title=''),fill='none')+
  scale_y_continuous(labels = scales::percent_format())+
  #facet_grid(vars(sub.group),scales="free_x",drop=TRUE,space="free")
  facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)
  
topptx(p, filename = "fig1.density_plot.pptx", width = 7, height = 3)

#保存为png
print(p)
graph2png(file=paste0("fig1.density_plot.png") ,width=7,height=5)

##Faceted plotting: Separately draw plots for the results from different tools-------------------------------
p=ggplot(input.data, aes(x=distance,weight=fren, 
                       fill=group,
                       color=group)) +
  geom_density(size=1,alpha=0.3)+scale_x_continuous(limits = c(-200, 200))+
  p_theme+
  scale_fill_manual(values=tool.col)+
  scale_color_manual(values=tool.col)+
  
  labs(x="Distance (bp)",y="Density")+
  theme(axis.text.x = element_text(angle = 0,hjust=0.6),
        strip.text.y = element_text(size = 8))+
  #guides(fill=guide_legend(title=NULL))+
  guides(fill='none',color="none")+
  scale_y_continuous(labels = scales::percent_format())+
  facet_grid(vars(group),scales="free_x",drop=TRUE,space="free")
topptx(p, filename = "fig1.sig_density_plot.pptx", width = 5, height = 5.5)

#保存为png
print(p)
graph2png(file=paste0("fig1.sig_density_plot.png") ,width=5,height=5.5)



##plot3.The proportion of overlapping pA sites-------------------------------------
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
pac.distance=readRDS("G:/scAPA/summary/temp/fig1/get_pac.distance(coord).rds")

sample.data$Sample_ID <- factor(sample.data$Sample_ID,
                                levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                         "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","GSM2333586","GSM2333587", "GSM2906396","GSM2906399",
                                         "3k","4k","5k","6k","GSM4712885","GSM4712907",
                                         "GSM3490689","GSM3490690","GSM3490691"),
                                labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                          "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                          "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                           "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                           "Arabidopsis1","Arabidopsis2","Arabidopsis3"))
pac.distance$srrID <- factor(pac.distance$srrID,
                             levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                      "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","GSM2333586","GSM2333587", "GSM2906396","GSM2906399",
                                      "3k","4k","5k","6k","GSM4712885","GSM4712907",
                                      "GSM3490689","GSM3490690","GSM3490691"),
                             labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                                        "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                        "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                                        "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                        "Arabidopsis1","Arabidopsis2","Arabidopsis3"))


samples=unique(sample.data$Sample_ID)

p.overlap.bar <- list()
pac.data <- data.frame()


#cutoff <- c(5)
#for(cutid in cutoff ){
  for(id in 1:length(samples)){
    
    sample_id <- samples[id]
    species <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Species)
    Protocol<- unique(sample.data[sample.data$Sample_ID==sample_id,]$Protocol)
    
    
    input.data <- subset(pac.distance,srrID== sample_id )
    input.data$type <- "Novel"
    input.data$type[which(  abs(input.data$distance) <=50)] <- "Overlapping"
    #input.data$label <- paste0(input.data$group,"_",input.data$type)
    input.data$label <- paste0(input.data$group,":",input.data$type)
    
    
    #计算每种工具识别的pA位点与已知pA重叠(距离<=50bp)和未重叠(novel)的数量
    pas.table <- as.data.frame(table(input.data$label))
    
    pas.table$group <- as.character(limma::strsplit2(pas.table$Var1,":")[,1])
    
    pas.table$type <- as.character(limma::strsplit2(pas.table$Var1,":")[,2])
    
    pas.table$species <- species
    pas.table$Protocol <- Protocol
    pas.table$srr <- sample_id
    pas.table <- pas.table%>%
      dplyr::group_by(group) %>% 
      dplyr::mutate(total=sum(Freq))
    
    pas.table$group <- factor(pas.table$group,
                              levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
    
    ## Plot PAC level
    p1.result <-plotNovel(novel.data=pas.table,input.level="pA Number")
    p1<-p1.result$p+labs(title=sample_id)
    #temp.label <- paste0("PAC",cutid)
    #p.overlap.bar[[temp.label]][[input.ID]] <- p1
    p.overlap.bar[[sample_id]] <- p1
    #p1.result$novel.data$level <-  temp.label
    pac.data <- plyr::rbind.fill(pac.data,p1.result$novel.data)
    
    cat(paste0(sample_id,": ",id,"/",length(samples)))
  }
#}

print(p.overlap.bar[[4]])

table(pac.data$type)

saveRDS(pac.data,file = "G:/scAPA/summary/temp/fig1/pac.data_utr3.refpa(coord).rds")

##（3.1）Plot the heatmap----------------------------------------------------------------
pac.data=readRDS("E:/scAPA/summary/temp/fig1/pac.data_utr3.refpa(coord).rds")



###在热图左侧同时添加样本和分组标签块--------------------------------------------

plot.Heatmap <- function(pac.data.need = NULL,input.label = "Real 3UTR pA%",
                         input.label2 =bquote("# 3UTR pA (x 10"^3~")"),
                         cols=NULL,
                         sample.info=NULL){
  
  p<- ggplot(pac.data.need, aes(y = factor(srr),
                                x = factor(group))) +        ## global aes
    geom_tile(aes(fill = rectheat),color="black",fill="#F5F5F5") + 
    #      ## to get the rect filled
    geom_point(aes(fill = pro.value, 
                   size =total),  shape = 21,colour= "#C7B19C",stroke = 1.2)  +    ## geom_point for circle illusion;shape=21:实心圆形；shape=23,实心菱形
    scale_fill_gradientn(
      colours = cols,
      name=input.label)+
    #scale_fill_gradientn(colours=c("#023FA5","#7D87B9","#BEC1D4",
    #                                  rev(colorspace::heat_hcl(20))),name=input.label)+       ## color of the corresponding aes
    #scale_size(range = c(1, 20))+             ## to tune the size of circles
    theme_bw()+
    theme(panel.background = element_rect(fill =  "gray100", color = 'white'),
          panel.grid.major = element_line(color = 'white', linetype = 'dotted'),
          panel.grid.minor = element_line(color = 'green', size = 2),
          text=element_text(family="Arial" ),
          axis.text.y =element_blank(),
          axis.title.x =element_text(color="black",size=9,family="Arial" ) ,
          #axis.title.y =element_text(color="black",size=12,family="Helvetica" ) ,
          axis.text.x =element_text(color="black",size=9,family="Arial" ,angle = 45,hjust=1) ,
          #axis.text.y =element_text(color=rep(c("black","blue","darkgreen"),c(9,5,2)),size=10,family="Helvetica" ) ,
          legend.text =element_text(color="black",size=9,family="Arial"),
          legend.title=element_text(color="black",size=10,family= "Arial" ))+
    #guides(size=guide_legend(title=input.label))+
    guides(size=guide_legend(title=input.label2))+
    labs(x=NULL,y=NULL)+
    #scale_size_continuous(range=c(4, 8))+
    geom_text(
      aes(label=round(pro.value,0)),
      size=3,
      #nudge_x=0, nudge_y=0.1,size=2.5,
      check_overlap=T,color="black"
    )
  
  temp <- pac.data.need[!duplicated(pac.data.need$srr),]
  temp <- temp[order(temp$srr),]
  print(temp$srr)
  temp$anno <- "Sample_ID"
  
  p1.sample <-ggplot(  temp, aes(y=srr, x=anno, fill=srr)) + geom_tile() + 
    scale_y_discrete(position="left") +
    scale_fill_manual(breaks = temp$srr,
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
  
  
  p1.special <- ggplot(temp, aes(y=srr, x=anno2, fill=sub.species)) + geom_tile() + 
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
  
  result <- list(p=p,p1=p1.special,p2=p1.sample)
  
  return(result)
  
}


####仅在热图左侧添加样本标签----------------------------------------------------
plot.Heatmap2 <- function(pac.data.need = NULL,input.label = "Real 3UTR pA%",
                         input.label2 =bquote("# 3UTR pA (x 10"^3~")"),
                         cols=NULL,
                         sample.info=NULL){
  
  p<- ggplot(pac.data.need, aes(y = factor(srr),
                                x = factor(group))) +        ## global aes
    geom_tile(aes(fill = rectheat),color="black",fill="#F5F5F5") + 
    #      ## to get the rect filled
    geom_point(aes(fill = pro.value, 
                   size =total),  shape = 21,colour= "#C7B19C",stroke = 1.2)  +    ## geom_point for circle illusion;shape=21:实心圆形；shape=23,实心菱形
    scale_fill_gradientn(
      colours = cols,
      name=input.label)+
    #scale_fill_gradientn(colours=c("#023FA5","#7D87B9","#BEC1D4",
    #                                  rev(colorspace::heat_hcl(20))),name=input.label)+       ## color of the corresponding aes
    #scale_size(range = c(1, 20))+             ## to tune the size of circles
    theme_bw()+
    theme(panel.background = element_rect(fill =  "gray100", color = 'white'),
          panel.grid.major = element_line(color = 'white', linetype = 'dotted'),
          panel.grid.minor = element_line(color = 'green', size = 2),
          text=element_text(family="Arial" ),
          axis.text.y =element_blank(),
          axis.title.x =element_text(color="black",size=9,family="Arial" ) ,
          #axis.title.y =element_text(color="black",size=12,family="Helvetica" ) ,
          axis.text.x =element_text(color="black",size=9,family="Arial" ,angle = 45,hjust=1) ,
          #axis.text.y =element_text(color=rep(c("black","blue","darkgreen"),c(9,5,2)),size=10,family="Helvetica" ) ,
          legend.text =element_text(color="black",size=9,family="Arial"),
          legend.title=element_text(color="black",size=10,family= "Arial" )
          #legend.position = "top"
          )+
    #guides(size=guide_legend(title=input.label))+
    guides(size=guide_legend(title=input.label2
    ))+
    labs(x=NULL,y=NULL)+
    #scale_size_continuous(range=c(4, 8))+
    geom_text(
      aes(label=round(pro.value,0)),
      size=3,
      #nudge_x=0, nudge_y=0.1,size=2.5,
      check_overlap=T,color="black"
    )
  
  temp <- pac.data.need[!duplicated(pac.data.need$srr),]
  temp <- temp[order(temp$srr),]
  print(temp$srr)
  temp$anno <- "Sample_ID"
  
  p1.sample <-ggplot(  temp, aes(y=srr, x=anno, fill=srr)) + geom_tile() + 
    scale_y_discrete(position="left") +
    scale_fill_manual(breaks = temp$srr,
                      values=species.col,
                      name="Sample_ID")+
    theme_minimal() + 
    theme(axis.text.x = element_blank(), 
          #axis.text.x = element_text(angle = 45,hjust=1),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank())+
    #axis.text.y =element_text(color="black",size=10,family="Helvetica" )) +
    xlab(NULL) + ylab(NULL) +guides(fill="none")
  
  
  result <- list(p=p,p1=p1.sample)
  
  return(result)
  
}


########################################################################




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
#pac.data.need <- subset(pac.data,type=="Overlapping"&level=="PAC5")
pac.data.need <- subset(pac.data,type=="Overlapping")
pac.data.need$pro.value <- round(pac.data.need$Freq/pac.data.need$total*100,2)
pac.data.need$rectheat <- 1
#pac.data.need$Protocol <- factor(pac.data.need$Protocol,
#                                 levels =sample.infor$Protocol)

##Calculate the average and tool order---------------
order=pac.data.need %>% 
  group_by( group) %>%  
  summarise(mean.pro = mean(pro.value, na.rm = TRUE)) %>%  
  ungroup() %>%
  arrange(desc(mean.pro))
#order.tool=order$group

#max_abs<- max(abs(res))



data=pac.data.need %>% 
  group_by(sub.species, group) %>%  
  summarise(pro.value = mean(pro.value, na.rm = TRUE),
            total = sum(total, na.rm = TRUE),
            Freq = sum(Freq, na.rm = TRUE)
            ) %>%  
  ungroup() %>%
  arrange(sub.species,desc(pro.value))

data$srr=data$sub.species
data$rectheat <- 1

data$group =factor(data$group,
                            levels =order$group )

#pac.data.need <- pac.data.need[order(pac.data.need$srr),]
#pac.data.need$total <- pac.data.need$total/10^3
data$total <- data$total/10^3


#col=colorRampPalette(c("#663d74","#fac03d"))(8)
p <- plot.Heatmap2(pac.data.need =data,
                  input.label = "Real 3UTR pA%",
                  input.label2 =bquote("# 3UTR pA (x 10"^3~")"),
                  cols=heatmap.col,
                  sample.info=NULL)

p$p %>% insert_left(p$p1,width=.1)
p.pac<-(p$p +scale_size_continuous(range=c(2,7)))  %>% insert_left(p$p1,width=.1) 

graph2png(file="fig1.utr3_pac_heatmap2(sub.species).png",
          width=6,height=4)
topptx(p.pac, filename = "fig1.utr3_pac_heatmap2(sub.species).pptx", width = 6, height = 4)


##（3.2）Calculate the proportion of overlapping pAs within the range of 10 to 100 bp----------------------------------
pac.distance$group <- factor(pac.distance$group,
                             levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))

cutoff.distance=seq(from = 10, to = 100, by = 10)
input.need=data.frame()
for(id in 1: length(cutoff.distance)){
  cutoff<- cutoff.distance[id]
  input.data <- pac.distance
  input.data  <- input.data[,c(1,2,4,5,6)]
  
  input.data$count <- 1
  
  input.data$type <- "Novel"
  input.data$type[which(  abs(input.data$distance) <=cutoff)] <- "Overlapping"
  
 
  input.data <- input.data %>%
    dplyr::group_by(species,protocol,group,type) %>% 
    dplyr::summarise(count=sum(count))
  
  #计算每组中pa的总数量
  input.data <- input.data %>%
    dplyr::group_by(species,protocol,group) %>% 
    dplyr::mutate(total.count=sum(count))
  
  input.data$fren <-   round(input.data$count/input.data$total.count*100,0)
  input.data$cutoff<-cutoff
  input.data<- as.data.frame(input.data)
  input.need<- rbind(input.need,input.data)
}

input.need=subset(input.need,input.need$type=="Overlapping")
input.need$sub.group=paste0(input.need$species,"_",input.need$protocol)
input.need$sub.group <- factor(input.need$sub.group,
                               levels =c("Human_10x","Mouse_10x","Arabidopsis_10x" ,"Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq"),
                               labels =c("Human (10X)","Mouse (10X)","Arabidopsis (10X)" ,"Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)") )


write.csv(input.need,file="G:/scAPA/summary/temp/fig1/tool.precision.csv",row.names = FALSE)

p=ggplot(input.need, aes(x = cutoff,y = fren, color = group,fill= group)) + 
  #geom_line(linetype="dashed",size=1) +
  geom_line() +
  geom_point(size=2,shape=21) +  
  scale_color_manual(values=tool.col)+
  scale_fill_manual(values=tool.col)+
  labs(x = "Cutoff(bp)", y = "Prediction Precision(%)")+p_theme+
  #guides(fill="none")+
  scale_x_continuous(limits = c(10, 100),breaks = seq(10,100,10))+
  theme(legend.title =element_blank(),
        #legend.position = c(.90, .3)
  axis.text.x = element_text(angle = 0,vjust=0.6),
  #strip =FALSE, 
  strip.text.x = element_text(size=9),
  element_rect(fill = "grey"))+
  guides(color=guide_legend(title=''),fill='none')+
  #facet_grid(vars(type),scales="free_x",drop=TRUE,space="free")
  facet_wrap(~ sub.group,scales="free_x",nrow=2,drop=TRUE)
  
topptx(p,filename ="fig1.pac_precision.pptx",width=7,height=5 )

print(p)
graph2png(file="fig1.pac_precision.png",width=7,height=5)


###plot4.genomic features of pAs and poly(A) signal--------------------------

##(1)load data----------------------
tool.list <- readRDS("G:/scAPA/summary/temp/fig1/gr_paclists_addDistance(coord).rds")
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)


annotation.lists<-readRDS("G:/scAPA/summary/temp/fig1/annotation.lists.rds")


#===Loading genomes

bsgenome.list <- list()
bsgenome.list[["Human"]] <- BSgenome.Hsapiens.UCSC.hg38
bsgenome.list[["Mouse"]] <- BSgenome.Mmusculus.UCSC.mm10
bsgenome.list[["Arabidopsis"]] <- BSgenome.Athaliana.TAIR.TAIR9

pas.annotation <- list()
data.stastic.all <- data.frame()
poly.signal.all <- data.frame()

samples=unique(sample.data$Sample_ID)

for(id in 1:length(samples)) {
  sample_id <- samples[id]
  species <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Species)
  
  Protocol <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Protocol)
  #temp <- tool.list[[srrID]]
  for(tools in names(tool.list[[sample_id]])){
    
    pac <- tool.list[[sample_id]][[tools]]$pac
   
    pac.result <- get.annotation(input.data = pac,
                                 annotation = annotation.lists[[ species ]],
                                 source=sample_id,
                                 cutoff.distance=50,
                                 species=species,
                                 ext3UTRlen=NULL,
                                 bsgenome =bsgenome.list[[species]])
    
    pas.annotation[[sample_id]][[tools]] <- pac.result$data.PACds
    
    #pac5.temp <- subset(pac.result$data.PACds@anno,count>=5)
    
    pac.temp=pac.result$data.PACds@anno
    
    ##统计基因组不同区域的pA数量
    temp <- data.stastic.temp(input.data=pac.temp)
    temp$source <- sample_id
    temp$species <- species
    temp$Protocol <-  Protocol
    temp$group <- tools
    #temp$pa_level <- 'PAC5'
    data.stastic.all <- plyr::rbind.fill(data.stastic.all,temp)
    
    #统计3UTR区域的poly(A)信号的分布
    utr3_pac.temp=subset(pac.result$data.PACds@anno,pac.result$data.PACds@anno$ftr=="3UTR")
    pa.label <- as.data.frame( table(paste0(utr3_pac.temp$pA.signal,
                                            "_",utr3_pac.temp$ftr,"_",
                                            utr3_pac.temp$type_distance) ))
    colnames(pa.label) <- c("Features1","No")
    
    pa.label$pa.signal <- as.character(limma::strsplit2(pa.label$Features1,"\\_")[,1])
    pa.label$Features <- as.character(limma::strsplit2(pa.label$Features1,"\\_")[,2])
    pa.label$type <- as.character(limma::strsplit2(pa.label$Features1,"\\_")[,3])
    pa.label$source <- sample_id
    pa.label$species <- species
    pa.label$Protocol <-  Protocol
    pa.label$group <- tools
    #pa.label$pa_level <- 'PAC5'
    poly.signal.all <- plyr::rbind.fill(poly.signal.all,pa.label)
  }
  
}

saveRDS(pas.annotation,file="G:/scAPA/summary/temp/fig1/gr_pac_addAnnotation(coord).rds")
saveRDS(data.stastic.all,file="G:/scAPA/summary/temp/fig1/genomicFeature_table.rds")
saveRDS(poly.signal.all,file="G:/scAPA/summary/temp/fig1/polyAsignal_table.rds")


##(2)Plotting genomic features===========================================
#

data.stastic.all<-readRDS("G:/scAPA/summary/temp/fig1/genomicFeature_table.rds")
sample.data<- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
#data.stastic.all$group <- factor(data.stastic.all$group,
#                                 levels=rev(c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
#                                            )


##(2.1) Plotting all pA at 3UTR, CDS,5UTR, intron

data.stastic.all.pa <-data.stastic.all %>% 
  dplyr::group_by(species,Protocol,source,group,Features) %>% 
  dplyr::summarise(Features.num=sum(No),.groups = 'keep') 


data.stastic.all.pa <-data.stastic.all.pa %>% 
  dplyr::group_by(species,Protocol,source,group) %>% 
  dplyr::mutate(pro_feature=Features.num/sum(Features.num)) 

data.stastic.all.pa$pro_3utr <-data.stastic.all.pa$pro_feature
data.stastic.all.pa$pro_3utr[-which(data.stastic.all.pa$Features=="3UTR")] <- NA


#===using mean value===============================
data.stastic.all.pa.average <- data.stastic.all.pa %>% 
  dplyr::group_by(species,group,Features) %>% 
  dplyr::summarise(pro_feature_mean=mean(pro_feature),.groups = 'keep') 

data.stastic.all.pa.average  %>% 
  dplyr::group_by(species,group) %>% 
  dplyr::summarise(test=sum(pro_feature_mean))

data.stastic.all.pa.average <- data.stastic.all.pa.average[order(data.stastic.all.pa.average$group,
                                                                 #data.stastic.all.pa.average$pa_level,
                                                                 data.stastic.all.pa.average$Features,
                                                                 decreasing = T),]

data.stastic.all.pa.average$sub.group=paste0(data.stastic.all.pa.average$species,"_",data.stastic.all.pa.average$group)


df_cumsum_mean <- ddply(data.stastic.all.pa.average, .(sub.group),
                        transform, 
                        label_ypos=cumsum(pro_feature_mean) - 0.5*pro_feature_mean)
df_cumsum_mean$label_ypos[-which(df_cumsum_mean$Features=="3UTR")] <- NA

df_cumsum_mean$species<- factor(df_cumsum_mean$species,
                                levels= c("Human","Mouse","Arabidopsis"))


##tool order
order<- data.stastic.all.pa[data.stastic.all.pa$Features=="3UTR",] %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(pro_feature_mean=mean(pro_feature),.groups = 'keep') %>%
  ungroup()%>%
  arrange(desc(pro_feature_mean))

df_cumsum_mean$group=factor(df_cumsum_mean$group,
                            levels = rev(order$group))


p=ggplot(data=df_cumsum_mean,
       aes(x=group, y=pro_feature_mean, fill=Features)) +
  geom_bar(stat="identity")+p_theme+
  labs(x=NULL,y="Percentage of genomic features (%)")+
  guides(fill=guide_legend(title=NULL))+
  #scale_fill_manual(values=wes_palette(n=5, name="Cavalcanti1"))+
  #scale_fill_manual(values=c("#299D8F","#E9C46A","#D87659","#b5182b","#934b43"))+
  scale_fill_manual(values=feature.col)+
  geom_text(aes(y=label_ypos, label=paste0(round(pro_feature_mean,2)*100,"%")), vjust=1, 
            color="white", size=3.5)+
  coord_flip()+
theme(legend.position = "top",
  axis.text.x = element_text(angle = 0,vjust=0.6),
  #strip =FALSE, 
  strip.text.x = element_text(size=9),
  element_rect(fill = "grey"))+
  #guides(color=guide_legend(title=''),fill='none')+
  scale_y_continuous(labels = scales::percent_format())+
  facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)

topptx(p, filename = "fig1.genomic_feature_all.mean.pptx", width = 7, height = 4)

print(p)
graph2png(file="fig1.genomic_feature_all.mean.png",width=7,height=4)


##(3)Plotting poly(A) signal===========================================
#
#
#
poly.signal.all <-readRDS("G:/scAPA/summary/temp/fig1/polyAsignal_table.rds")
#poly.signal.all$group <- factor(poly.signal.all$group,
#                                levels=rev(c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
#                                           )



##(3.1) All signals
#不同样本数据中，各工具识别的3UTR pA附近不同polya信号的数量
poly.signal.all.pa <- poly.signal.all  %>% 
  dplyr::group_by(species,Protocol,source,group,pa.signal) %>% 
  dplyr::summarise(Features.num=sum(No),.groups = 'keep') 

#不同样本数据中，各工具识别的pA附近不同polya信号的比例
poly.signal.all.pa  <- poly.signal.all.pa %>% 
  dplyr::group_by(species,Protocol,source,group) %>% 
  dplyr::mutate(pro_feature=Features.num/sum(Features.num)) 

poly.signal.all.pa$pa.signal <- factor(poly.signal.all.pa$pa.signal,
                                       levels=rev(c("AATAAA","1Variants","Others")) ,
                                       labels=rev(c("AATAAA","1nt-variants","No signals")))


#=============getting average==================================================

poly.signal.all.pa.average  <- poly.signal.all.pa %>% 
  dplyr::group_by(species,group,pa.signal) %>% 
  dplyr::summarise(pro_feature_mean=mean(pro_feature),.groups = "keep") 


poly.signal.all.pa.average <- poly.signal.all.pa.average[order(poly.signal.all.pa.average$group,
                                                               #poly.signal.all.pa.average$pa_level,
                                                               poly.signal.all.pa.average$pa.signal,
                                                               decreasing = T),]

poly.signal.all.pa.average$sub.group=paste0(poly.signal.all.pa.average$species,"_",poly.signal.all.pa.average$group)


df_cumsum_mean <- ddply(poly.signal.all.pa.average, .(sub.group),
                        transform, 
                        label_ypos=cumsum(pro_feature_mean) - 0.5*pro_feature_mean)


df_cumsum_mean$species<- factor(df_cumsum_mean$species,
                                levels= c("Human","Mouse","Arabidopsis"))


##tool order
order<- poly.signal.all.pa[poly.signal.all.pa$pa.signal=="No signals" & !poly.signal.all.pa$species=="Arabidopsis",] %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(pro_feature_mean=mean(pro_feature),.groups = 'keep') %>%
  ungroup()%>%
  arrange(desc(pro_feature_mean))

df_cumsum_mean$group=factor(df_cumsum_mean$group,
                            levels = order$group)


#class(df_cumsum_mean$pa.signal)
p=ggplot(data=df_cumsum_mean,
       aes(x=group, y=pro_feature_mean, fill=pa.signal)) +
  geom_bar(stat="identity")+p_theme+
  labs(x=NULL,y="Percentage of poly(A) signals (%)")+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=rev(c("#992813","#a24f47","#eccab7")))+
  #scale_fill_manual(values=rev(c("#F2AD00","#ECCBAE","#D6D6D6")))+
  geom_text(aes(y=label_ypos, label=paste0(round(pro_feature_mean,2)*100,"%")), vjust=1, 
            color="white", size=2.5)+
  coord_flip()+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 0,vjust=0.6),
        #strip =FALSE, 
        strip.text.x = element_text(size=9),
        element_rect(fill = "grey"))+
  guides(fill = guide_legend(reverse=TRUE,title=""))+
  scale_y_continuous(labels = scales::percent_format())+
  facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)


topptx(p, filename = "fig1.utr3_polyA_signal_all.mean.pptx", width = 7, height = 4)

print(p)
graph2png(file="fig1.utr3_polyA_signal_all.mean.png",width=7,height=4)

##plot5.Single nucleotide distribution profile---------------------------------------------------------

#load gedata
bsgenome.list <- list()
bsgenome.list[["Human"]] <- BSgenome.Hsapiens.UCSC.hg38
bsgenome.list[["Mouse"]] <- BSgenome.Mmusculus.UCSC.mm10
bsgenome.list[["Arabidopsis"]] <- BSgenome.Athaliana.TAIR.TAIR9

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

temp.path="G:/scAPA/summary/temp/fig1/fafile"
setwd(temp.path)

GSElists=unique(sample.data$Sample)

fafiles.utr3 <- c()
for(s in 1: length(GSElists)){
  
  sample=GSElists[s]
  species <- unique(sample.data[sample.data$Sample==sample,]$Species)
  
  protocol <- unique(sample.data[sample.data$Sample==sample,]$Protocol)
  
  GSMlists=unique(sample.data[sample.data$Sample==sample,]$Sample_ID)
  
  all.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  
  
  
  for(i in 1:length( GSMlists)) {

  sample_id=GSMlists[i]
  tool.lists <- names(all.pac[[sample_id]])
  
  for(j in 1:length(tool.lists)){
    tools=tool.lists[j]
    bsgenome=bsgenome.list[[species]]
    pac <-all.pac[[sample_id]][[tools]]
    
    if(species=="Arabidopsis"){
    
     
      pac@anno$chr=paste0("Chr",pac@anno$chr)
    }
    
    fafiles.all.pac <- faFromPACds(PACds= pac,
                                    bsgenome=  bsgenome,
                                    what = "updn", fapre = paste0(sample_id,"_",tools,"_200nt"),
                                    up = -100, dn = 100, byGrp ="ftr", chrCheck = FALSE)
    
   
    ###keep 3utr sequences
    fafiles.temp <-paste0(sample_id,"_",tools,"_200nt.3UTR.fa")
    fafiles.utr3  <- c(fafiles.utr3, fafiles.temp)
    
    message(paste0(sample_id,": ",j,"/",length(tool.lists)))
  }
  }
}



fafiles.utr3
  

#select.fa <- fafiles.utr3[grep(paste0("MI_DAY3_TIP+SHAM_TIP"),fafiles.utr3)]
#select.fa <- fafiles.utr3[1:8]
#when sample data is Arabidopsis
#select.fa <- fafiles.utr3[1:10]

samples=unique(sample.data$Sample_ID)
labels=c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
         "Mouse_Bone1","Mouse_Bone2","Mouse_Hypothalamus1", "Mouse_Hypothalamus2",
         "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
         "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
         "Arabidopsis1","Arabidopsis2","Arabidopsis3")

temp.path="G:/scAPA/summary/temp/fig1/fafile"
setwd(temp.path)

for(i in 1: length(samples)){
  sample=samples[i]
  
  tool.lists=unique(sample.data[sample.data$Sample_ID==sample,]$Tool)
  
  fa=c()
  for(j in 1:length(tool.lists)){
    tool=tool.lists[j]
    temp.fa=paste0(sample,"_",tool,"_200nt.3UTR.fa")
    fa=c(fa, temp.fa)
    
  }
  
  p<-plotATCGforFAfile.cp(faFiles=fa, ofreq = TRUE, opdf = F, refPos = 101,
                          filepre = "",mergePlots = T)
  print(p$ATCGpicture)
  
  din <- p$din   
  din$source.raw <- din$source
  
  if(sample=="MI_DAY3_TIP+SHAM_TIP"){
    
    din$source  <- as.character(limma::strsplit2(din$source.raw,"\\_")[,5])
    
    din$number <- as.integer(limma::strsplit2(din$source.raw,"\\#")[,2])
    
  }else if(sample=="MI_DAY7_GFP" | sample=="SHAM_DAY7_GFP"){
  
    din$source  <- as.character(limma::strsplit2(din$source.raw,"\\_")[,4])
    
    din$number <- as.integer(limma::strsplit2(din$source.raw,"\\#")[,2]) 
    
  } else{
    
    din$source  <- as.character(limma::strsplit2(din$source.raw,"\\_")[,2])
    
    din$number <- as.integer(limma::strsplit2(din$source.raw,"\\#")[,2]) 
    
  }
  
  
  din$source <- factor(din$source,
                       levels=c("MAAPER","polyApipe","scUTRquant","scAPAtrap","Sierra","scraps","SCAPTURE","Infernape","SCAPE","scAPA"),
                       labels = c("MAAPER","polyApipe","scUTRquant","scAPAtrap","Sierra","scraps","SCAPTURE","Infernape","SCAPE","scAPA"))
  
  
  
  #output <- "G:/scAPA/summary/Figure/fig1"
  #if(!dir.exists(output)){
  #  dir.create(output)
  #}
  #setwd(output)
  
  #sample.title="Arabidopsis2"
  sample.title=labels[i]
  
  p2=ggplot(din, aes(x=pos, y=freq, group=base)) +
    geom_line(aes(colour=base)) +
    labs(x="Position",y='Base component',title= sample.title) + p_theme+
    #facet_wrap(~source,scales="free_x",nrow=2,drop=TRUE)+
    facet_wrap(~source,scales="free_x",nrow=2,drop=TRUE)+
    scale_color_manual(values = base.col)+ 
    #scale_color_brewer(palette="Set1")+
    guides(color=guide_legend(title=NULL))+
    theme(#legend.position  ="none",
      axis.text.x = element_text(angle = 0,vjust=0.6,size=8),
      #strip =FALSE,
      strip.text.x = element_text(size=8))
  
  topptx(p2, filename =paste0("G:/scAPA/summary/Figure/fig1/nucleotide.fre/nucleotide frequency(",sample.title,").pptx"), width = 7.5, height = 4)
  
  print(p2)
  graph2png(file=paste0("G:/scAPA/summary/Figure/fig1/nucleotide.fre/nucleotide frequency(",sample.title,").png"),width=7.5,height=4)
  
  
  print(paste0("Finish plot :",sample," (",i,"/",length(samples),")"))
  
}



##plot6.Chi-square test----------------------------------------------------------------
 

####select 5000 3UTR PAC from reference known PAC as gold standard
#din.object.ref <- readRDS("G:/R/code/scAPA_Demo_figure/Ref_select_50003UTR.Rdata")
#din.object <- din.object.ref$din.object
setwd("G:/scAPA/summary/Figure/fig1")
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

## input fasta file generated from  mouse_sperm data
#fasta.file <- list.files("G:/scAPA/summary/temp/fig1/fafile","GSM2803334\\S+_200nt.3UTR.fa")
#fasta.file <- list.files("G:/scAPA/summary/temp/fig1/fafile","\\S+_200nt.3UTR.fa")
#fasta.file 

####(1)Calculate the distribution frequency of bases----------------------------------------------------------

GSMlists=unique(sample.data$Sample_ID)
seq.num=1500

din.object.list <- list()

for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  sample_id<-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==sample_id,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==sample_id,]$Sample)
 
  #sample_id <- as.character(limma::strsplit2(fafile,"\\_")[,1])
  #tool <- as.character(limma::strsplit2(fafile,"\\_")[,2])
  #pa.level <- as.character(limma::strsplit2(fafile,"\\_")[,3])
  
  tool.lists=sample.data[sample.data$Sample_ID==sample_id,]$Tool
  
  for(j in 1:length(tool.lists)){
    tool=tool.lists[j]
  
    fafile=paste0(sample_id,"_",tool,"_200nt.3UTR.fa")
 
  seq=readDNAStringSet(paste0("G:/scAPA/summary/temp/fig1/fafile/",fafile), format="fasta")
  nseq=length(seq)
  if(nseq<seq.num){
    next
  }else{
   
    din<-consensusMatrix(seq, as.prob=TRUE, baseOnly=TRUE)
    din<- as.data.frame(t(din))
    din$other <- NULL
    din$pos <- c(1:nrow(din))
    din$number <- nseq
    din.data2 <- reshape2::melt(din[,1:5],id.vars='pos',variable.name = 'base', value.name = 'freq')
    ATCGpicture.all <-ggplot(din.data2,aes(pos,y=freq,colour=base)) + geom_line() +
      xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw() +
      #scale_color_manual(values = base.col)+
      ggtitle("All-Whole-Gene-Body")
    
    din.object.list[["allseq"]][[sample_id]][[tool]] <- din
    din.object.list[["allfig"]][[sample_id]][[tool]] <- ATCGpicture.all
    
    for(i in 1:100){
     
      index <- sample(1:nseq,size=seq.num,replace = F)
      seqsub <- seq[index]
      dinsub<-consensusMatrix(seqsub, as.prob=TRUE, baseOnly=TRUE)
      dinsub<- as.data.frame(t(dinsub))
      dinsub$other <- NULL
      dinsub$pos <- c(1:nrow(dinsub))
      dinsub$number <- length(seqsub)
      
      din.subdata <- reshape2::melt(dinsub[,1:5],id.vars='pos',variable.name = 'base', value.name = 'freq')
      ATCGpicture.sub <-ggplot(din.subdata,aes(pos,y=freq,colour=base)) + geom_line() +
        xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw() +
        #scale_color_manual(values = base.col)+
        ggtitle("1500-Whole-Gene-Body")
      
      din.object.list[["seq"]][[sample_id]][[tool]][[i]] <- dinsub
      din.object.list[["fig"]][[sample_id]][[tool]][[i]] <- ATCGpicture.sub
    }
  }
  message(paste0(sample_id,": ",j,"/",length(tool.lists))) 
  
  }
}

names(din.object.list[["fig"]])
names(din.object.list[["fig"]][["GSM2803334+GSM2803335"]])


saveRDS(din.object.list,file=paste0("G:/scAPA/summary/temp/fig1/",seq.num,"seq.din.obj.lists.rds"))

#=============================================================================
####(2) get chisq.test for collected pA---------------------------------

din.object=readRDS("E:/scAPA/summary/refPA/refseq.1500.3UTR.rds")
#din.object.ref <- readRDS("G:/R/code/scAPA_Demo_figure/Ref_select_50003UTR.Rdata")
#din.object <- din.object.ref$din.object

#din.object.list=readRDS(paste0("E:/scAPA/summary/temp/fig1/1500seq.din.obj.lists.rds"))

GSMlists=unique(sample.data$Sample_ID)
chiseq.object.list <-list()
res=data.frame()

for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  GSM <-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  ref.data <- din.object[[species]][["ref"]][["repseq"]]
  ref.data<- ref.data[,1:4]*ref.data$number
  ref.input <- c(ref.data[,1],
                 ref.data[,2],
                 ref.data[,3],
                 ref.data[,4])
  
  
      tool.lists<-names(din.object.list[["seq"]][[GSM]])
      
      for(j in 1 : length(tool.lists)) {
        tool =tool.lists[j]
        temp <- din.object.list[["seq"]][[GSM]][[tool]]
        
        for(n in 1:100){
          obj.data <- temp[[n]]
          obj.data<- obj.data[,1:4]*obj.data$number
          obj.input <- c(obj.data[,1],
                         obj.data[,2],
                         obj.data[,3],
                         obj.data[,4])
        
          x <- matrix(c(as.integer(ref.input) ,
                        as.integer(obj.input)),nrow=2,byrow=TRUE)  
          
          chiseq.object.list[[GSM]][[tool]][["statistic"]][[n]] <- chisq.test(x)$statistic
          chiseq.object.list[[GSM]][[tool]][["pvalue"]][[n]] <- chisq.test(x)$p.value
          
        
        }
        
        
         
        statistic <- as.numeric(unlist(chiseq.object.list[[GSM]][[tool]][["statistic"]]))
        pvalue <- as.numeric(unlist(chiseq.object.list[[GSM]][[tool]][["pvalue"]]))
        
        temp.res <- data.frame(species=species,
                           protocol=protocol,
                           sample=sample,
                           sample_id=GSM,
                           tool=tool,
                           statistic=statistic,
                           pvalue=pvalue)
        temp.res$id <- 1:100
        
        res=rbind(res,temp.res)
        
        message(paste0(sample,": ",j,"/",length(tool.lists))) 
      }
  }


res$label.value=res$statistic/201  
write.csv(res,file="E:/scAPA/summary/temp/fig1/1500seq.chisp.res.csv",row.names = FALSE)


#####(3)plot--------------------------------------------------------------
###Plot specified samples
res=read.csv("E:/scAPA/summary/temp/fig1/1500seq.chisp.res.csv")
res$sample_id <- factor(res$sample_id,
                        levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                                 "GSM2906396","GSM2906399","GSM2333586","GSM2333587","GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287",
                                 "3k","4k","5k","6k","GSM4712885","GSM4712907",
                                 "GSM3490689","GSM3490690","GSM3490691"),
                        labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2","Mouse_Bone1","Mouse_Bone2",
                                   "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                                   "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                                   "Arabidopsis1","Arabidopsis2","Arabidopsis3"))

res$tool <- factor(res$tool,
                   levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"),
                   labels = c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))
table(res$tool)

######（3.1）representative samples from different sequencing protocols-------------------------------
samples=c("Mouse_Sperm","Mouse_Bone1","Mouse_Hypothalamus1","Mouse_Intestine1","Human_Pbmc3","Arabidopsis3")

#for(i in unique(data.obj.fig$sample_id)){
  sub.res=subset(res,res$sample_id%in% samples)
  
  table(sub.res$sample_id)
  
  
  
  
  sub.res$sample_id <- factor(sub.res$sample_id,
                              levels=c("Mouse_Sperm","Human_Pbmc3","Arabidopsis3","Mouse_Intestine1","Mouse_Hypothalamus1","Mouse_Bone1"),
                              labels=c("Mouse_Sperm(10X)","Human_Pbmc3(10x)","Arabidopsis3(10x)",
                                       "Mouse_Intestine1(CEL-seq)","Mouse_Hypothalamus1(Drop-seq)","Mouse_Bone1(Microwell-seq)"))
  
  ###################################################################################################
  
  #####（3.2）different sequencing protocols---------------------------------------------------
  tech="Drop-seq"
  sub.res=subset(res,res$protocol %in% tech)
  table(sub.res$sample_id)
  
  #######################################################################################
  
  ###（3.3）Grouped by species---------------------------------------------------------
  res=read.csv("E:/scAPA/summary/temp/fig1/1500seq.chisp.res.csv")
  sub.res=res
  
  
  sub.res.mean=sub.res %>% 
    dplyr::group_by(species,tool) %>% 
    dplyr::summarise(statistic.mean=mean(statistic),.groups = "keep") %>%
    ungroup()%>%
    arrange(species,statistic.mean)
  
  
  #基于所有数据中卡方检验的平均值给工具排序
  order= sub.res %>%
    dplyr::group_by(tool) %>% 
    dplyr::summarise(statistic.mean=mean(statistic),.groups = "keep") %>%
    ungroup()%>%
    arrange(statistic.mean)
  
  order2= sub.res[!sub.res$species=="Arabidopsis",] %>%
    dplyr::group_by(tool) %>% 
    dplyr::summarise(statistic.mean=mean(statistic),.groups = "keep") %>%
    ungroup()%>%
    arrange(statistic.mean)
  
  sub.res$species=factor(sub.res$species,
                         level=c("Human","Mouse","Arabidopsis")
  )
  
  sub.res$tool=factor(sub.res$tool,
                      levels = rev(order$tool)
  )
  
  
  p=ggplot(sub.res, aes(x=tool,y=statistic/201,color=tool)) +
    geom_boxplot(width=0.7,outlier.shape = 23,outlier.size = 1,outliers=F)+
    scale_color_manual(values=tool.col)+
    stat_boxplot(geom='errorbar',width=0.15)+
    xlab(NULL)+labs(y="Chi-squared test statistic")+
    p_theme+
    guides(fill="none",color="none")+
    theme(
      strip.text.x = element_text(size=8),
      element_rect(fill = "grey"))+
         coord_flip()+
    #facet_wrap(~ species,scales="free_x",nrow=1,drop=TRUE)
    facet_wrap(~ species,nrow=1,drop=TRUE)
  
  topptx(p,filename =paste0("fig1.1500seq.chisquard.pptx"),width=7,height=5 )
  #topptx(p,filename =paste0("fig1.1500seq.chisquard_",tech,".pptx"),width=7,height=5 )
  
  print(p)
  graph2png(file=paste0("fig1.1500seq.chisquard.png"),width=7,height=5)
  
  


###plot7.single nucleotide profile of reference pA --------------------------------------------

  bsgenome.list <- list()
  bsgenome.list[["Human"]] <- BSgenome.Hsapiens.UCSC.hg38
  bsgenome.list[["Mouse"]] <- BSgenome.Mmusculus.UCSC.mm10
  bsgenome.list[["Arabidopsis"]] <- BSgenome.Athaliana.TAIR.TAIR9
  
  samples=c("Mouse","Human","Arabidopsis")
  
  refpa=readRDS("E:/scAPA/summary/temp/fig1/anno.utr3.ref.pacds.rds")
 
  setwd("E:/scAPA/summary/refPA/fa")
  
  fa=c()
  for(i in 1: length(samples)){
    sample=samples[i]
    genome=samples[i]
    
    PACds=refpa[[sample]]
    bsgenome=bsgenome.list[[sample]]
    
    if(sample=="Arabidopsis"){
  
      PACds@anno$chr=paste0("Chr",PACds@anno$chr)
    }
    
    temp.fa = faFromPACds(PACds, bsgenome
                          , what = "updn", fapre = paste0("refPA_",genome,"_200nt"),
                          up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
    
    fa=c(fa, temp.fa)
    
  }
  p<-plotATCGforFAfile.cp(faFiles=fa, ofreq = TRUE, opdf = F, refPos = 101,
                          filepre = "",mergePlots = T)
  print(p$ATCGpicture)
  
  
  din <- p$din   
  din$source.raw <- din$source
  
 
  din$source  <- as.character(limma::strsplit2(din$source.raw,"\\_")[,2])
  
  din$number <- as.integer(limma::strsplit2(din$source.raw,"\\#")[,2]) 
  
  din$source <- factor(din$source,
                       level=c("Human","Mouse","Arabidopsis"))
  
  sample.title=""
  
  p2=ggplot(din, aes(x=pos, y=freq, group=base)) +
    geom_line(aes(colour=base)) +
    labs(x="Position",y='Base component',title= sample.title) + p_theme+
    #facet_wrap(~source,scales="free_x",nrow=2,drop=TRUE)+
    facet_wrap(~source,scales="free_x",nrow=1,drop=TRUE)+
    scale_color_manual(values = base.col)+ 
    guides(color=guide_legend(title=NULL))+
    theme(#legend.position  ="none",
      axis.text.x = element_text(angle = 0,vjust=0.6,size=8),
      #strip =FALSE,
      strip.text.x = element_text(size=8))
  
  topptx(p2, filename =paste0("E:/scAPA/summary/Figure/fig1/nucleotide.fre/nucleotide frequency(refPA).pptx"),
         width = 8, height = 4)
  
  print(p2)
  graph2png(file=paste0("E:/scAPA/summary/Figure/fig1/nucleotide.fre/nucleotide frequency(refPA).png"),
            width=8,height=4)
  
  
  
