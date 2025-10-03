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
#library(BSgenome.Mmusculus.UCSC.mm10)

heatmap.col=c("#313772","#2c4ca0","#326db6","#478ecc","#75b5dc","#fee3ce","#eabaa1","#dc917b","#d16d5b","#c44438","#b7282e")


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

tool.col=c("scAPAtrap"="#064D4B","Sierra"="#fbd26a","SCAPE"="#ABC8E5","scAPA"="#6666CC","polyApipe"="#e8b9c5","MAAPER"="#64d9b8",
           "SCAPTURE"="#b5182b","Infernape"="#775360","scUTRquant"="#7fb80e","scraps"="#FF6B00")


tool.col2=c("#064D4B","#fbd26a","#b5182b","#775360","#ABC8E5","#6666CC","#64d9b8","#e8b9c5","#7fb80e","#FF6B00")


tool.col3=c("#064D4B","#fbd26a","#002fa7","#6666CC","#b56387","#64d9b8","#b5182b","#775360","#7fb80e","#FF6B00")

seq.col=c("10x"="#a5cc26","CEL-seq"="#A46ACD","Drop-seq"="#E78D35","Microwell-seq"="#002fa7")


##########Create a list of cell type annotation information------------------------------------------------------
GSElists=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
           "arab")

all_datalists=c("merge","merge","merge","merge","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399","GSM2333586","GSM2333587",
                "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
                "GSM3490689","GSM3490690","GSM3490691")

all_GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
               "GSM2333586","GSM2333587", "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490690","GSM3490691")


anno.lists=list()

for(s in 1: 11){
  
  sample=GSElists[s]
  GSE=paste0("G:/scAPA/barcode/",sample)
  
  
  if(s %in% 1:4){
    
    datalists=c("merge")
    GSMlists=all_GSMlists[s]
    ##（mouse_sperm;mouse_tcell;mouse_tip:mouse_esc)的GSM编号均为merge
    #GSMlists=dir(file.path)
  }else if (s ==8){
    GSMlists=c("GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287")
    datalists=GSMlists
  }else{
    GSMlists=dir(GSE)
    datalists=GSMlists
  }
  
  
  for(i in 1: length(GSMlists)){
    GSM=GSMlists[i]
    data=datalists[i]
    
    if(s==8){
      file.path=paste0(GSE)
    }else{
      file.path=paste0(GSE,"/",data)
      
    }
    
    temp.anno=read.csv(paste0(file.path,"/cell.csv"))
    
    anno.lists[[GSM]]=temp.anno
    
    message(paste0(GSM,": ",i,"/",length(GSMlists)))
  }
 } 
    

saveRDS(anno.lists,file="G:/scAPA/summary/temp/fig2/bech_pacount/anno.lists.rds") 


#############一：gene level:Correlation with gene count-------------------------------------------------------------------



GSElists=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
           "arab")

all_datalists=c("merge","merge","merge","merge","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399","GSM2333586","GSM2333587",
                "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
                "GSM3490689","GSM3490690","GSM3490691")


all_GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
               "GSM2333586","GSM2333587", "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490690","GSM3490691")


toolfile.lists=c("scapatrap","sierra","scape","scapa","polyapipe","maaper","scapture","infernape","scutrquant","scraps")
tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps")


#Load the gene expression matrix into the Seurat object------------------------------------------------------------------------------------

##Load the genome annotation file to enable conversion between gene symbol and gene emsembl

reference.file="G:/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf"

ref.gene<-parseGenomeAnnotation(reference.file)

gene.info=ref.gene[["anno.frame"]][,c("transcript_id","gene_id","gene_name")] 


#toolfile=toolfile.list[1]
#tool=tool.list[1]
#for(s in 1: length(GSElists)){
gcounts=list()

for(s in 1: 11){
  
  sample=GSElists[s]
  GSE=paste0("G:/scAPA/gene/",sample)
  
  
  if(s %in% 1:4){
   
    datalists=c("merge")
    GSMlists=all_GSMlists[s]
    ##四份合并的数据集（mouse_sperm;mouse_tcell;mouse_tip:mouse_esc)的GSM编号均为merge
    #GSMlists=dir(file.path)
  }else{
    GSMlists=dir(GSE)
    datalists=GSMlists
  }
  
  for(i in 1: length(GSMlists)){
    GSM=GSMlists[i]
    data=datalists[i]
    
    file.path=paste0(GSE,"/",data)
    
    if(data %in% c("GSM2906396","GSM2906399") ){
      count_symbol=readRDS(paste0(file.path,"/","gcount.rds"))
      info=intersect(row.names(count_symbol),gene.info$gene_name)
      count_symbol=count_symbol[info,]
      
      loc=match(row.names(count_symbol),gene.info$gene_name)
      table(is.na(loc))
      count_ensembl=count_symbol
      row.names(count_ensembl)=gene.info$gene_id[loc]
      
      count_symbol=as.matrix(count_symbol)
      count_ensembl=as.matrix(count_ensembl)
      
      saveRDS(count_symbol,file=paste0(file.path,"/","gcount.rds"))
      saveRDS(count_ensembl,file=paste0(file.path,"/","gcount_ensembl.rds"))
      
    }else{
      count_symbol=Seurat::Read10X(
        file.path,
        gene.column = 2,  
        strip.suffix = FALSE
      )
      
      count_ensembl=Seurat::Read10X(
        file.path,
        gene.column = 1, 
        strip.suffix = FALSE
      )
      
      saveRDS(count_symbol,file=paste0(file.path,"/","gcount.rds"))
      saveRDS(count_ensembl,file=paste0(file.path,"/","gcount_ensembl.rds"))
    }
    
    gcounts[[GSM]]=count_ensembl
    
    message(paste0(GSM,": ",i,"/",length(GSMlists),sample))
  }
}

saveRDS(gcounts,file="G:/scAPA/summary/gene/gcounts.rds") 

####1.Correlation between gene expression and the sum of APA transcript expression on the gene---------

#gcounts=readRDS("G:/scAPA/summary/gene/gcounts.rds")


cor.test="pearson" 
res=data.frame()
count.lists=list()
for(s in 1: length(GSElists)){
  
  sample=GSElists[s]
  
  species <- ifelse(s %in% 1:8, "Mouse",
                    ifelse(s %in% 9:10, "Human",
                           ifelse(s == 11, "Arabidopsis",NA)))
  
  protocol<- ifelse(s %in% c(1:5,9,10,11), "10x",
                    ifelse(s ==6, "Microwell-seq",
                           ifelse(s == 7, "Drop-seq",
                                  ifelse(s == 8, "CEL-seq", NA))))
  
  
  
  if(sample=="mouse_bone"){
    tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scraps")
    
    
  } else if(sample=="mouse_int"){
    tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","MAAPER","SCAPTURE","Infernape","scraps")
  }else if(sample=="arab"){
    tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","scraps")
    
  }else{
    tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps")
    
  }
  
  
  
  all.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  GSMlists=names(all.pac)
  
  
  for(i in 1: length(GSMlists)){
    GSM=GSMlists[i]
    
    gcount=gcounts[[GSM]]
    
    colnames(gcount)=gsub("-.*","",colnames(gcount))
    
    for(j in 1: length(tool.lists)){
      tool=tool.lists[j]
      
      ##Data that did not run successfully
      conditions <- list(
        c(GSM == "GSM3490689", tool == "SCAPE"),
        c(GSM == "GSM2906399", tool == "scAPA"),
        c(GSM == "3k", tool == "scUTRquant"),
        c(GSM == "6k", tool == "scUTRquant")
      )
      
      if(any(sapply(conditions, all))){
        print("PAC data not exist")
        next 
      }else{
        pac=all.pac[[GSM]][[tool]]
        
        if(sample=="arab" & tool=="SCAPE"){
          colnames(pac@counts)= gsub(paste0(GSM,"."),"", colnames(pac@counts))
          row.names(pac@colData)= gsub(paste0(GSM,"."),"",  row.names(pac@colData))
          pac@colData$barcode=row.names(pac@colData)
        }
        
        
        colnames(pac@counts)=gsub("-.*","",colnames(pac@counts))
        row.names(pac@colData)=gsub("-.*","",row.names(pac@colData))
        pac@colData$barcode=row.names(pac@colData)
        
        ##MAAPER only outputs total expression levels
        if(tool=="MAAPER"){
         
          info.bc=intersect(colnames(all.pac[[GSM]][["scAPAtrap"]]@counts),colnames(gcount))
          info.gene=intersect(pac@anno$gene,row.names(gcount))
          
          sub.gcount=gcount[info.gene,info.bc]
          
          
          #pac@counts=pac@counts[,info.bc]
          pac@anno$count=pac@counts[,1]
          pac@anno$pa_id=paste0(pac@anno$chr,":",pac@anno$strand,":",pac@anno$coord)
          pac.count=pac@anno[,c("pa_id","gene","count")]
          
          pac.count=subset(pac.count,pac.count$gene %in% info.gene)
          
          
        }else{
          
          #Keep barcodes and genes common to both the gene expression matrix 
          #and the APA transcript expression matrix
          info.bc=intersect(colnames(pac@counts),colnames(gcount))
          info.gene=intersect(pac@anno$gene,row.names(gcount))
          
          sub.gcount=gcount[info.gene,info.bc]
          
          ##Calculate the total expression levels of genes and all APA transcripts on those genes
          pac@counts=pac@counts[,info.bc]
          pac@anno$count=Matrix::rowSums(pac@counts)
          pac@anno$pa_id=paste0(pac@anno$chr,":",pac@anno$strand,":",pac@anno$coord)
          pac.count=pac@anno[,c("pa_id","gene","count")]
          
          pac.count=subset(pac.count,pac.count$gene %in% info.gene)
          
        }
        
        
        pac.res= pac.count %>% 
          dplyr::group_by(gene) %>% 
          dplyr::summarise(counts=sum(count))
        
        gene.res=data.frame(row.names = row.names(sub.gcount),
                            gene=row.names(sub.gcount),
                            counts=Matrix::rowSums(sub.gcount))
        
        gene.res=gene.res[pac.res$gene,]
        
        cor=psych::corr.test(pac.res$counts,gene.res$counts,  
                             use = "pairwise", 
                             method=cor.test, 
                             adjust = "fdr"   
        )
        
       count.temp=data_frame(row.names=pac.res$gene,
                             apa.counts=pac.res$counts,
                             gene.counts=gene.res$counts)
       count.lists[[GSM]][[tool]]=count.temp
        
        temp=data.frame(species=species,
                        sample=GSM,
                        protocol=protocol,
                        method=tool,
                        cor=cor$r,
                        p.adj=cor$p.adj,
                        gene.num=cor$n,
                        cell.num=length(info.bc),
                        cor.method=cor.test
        )
        
        res=rbind(res,temp)
      }
      
      message(paste0(sample,"(",GSM,")",": ",j,"/",length(tool.lists)))
    }                
    
  }
}

write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/gcount.cor.csv",row.names = FALSE)                   
                
saveRDS(count.lists,file ="G:/scAPA/summary/temp/fig2/bech_pacount/counts.lists.rds" )


#######plot----------------------------------------------------------------------

res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/gcount.cor.csv")
res$group=paste0(res$species,"_",res$protocol)

##(1)Grouped by different sequencing protocols-------------------------------------------
data=res %>%
  dplyr::group_by(protocol,method) %>%
  dplyr::summarise(
    mean_cor = mean(cor, na.rm = TRUE),
    #sd_value = sd(value, na.rm = TRUE),
    .groups = "keep"  
  )%>%
  arrange(protocol,method)

min(data$mean_cor)


data$protocol<-factor(data$protocol,
                   level=c("10x" ,"CEL-seq","Drop-seq","Microwell-seq"),
                   labels=c("10X" ,"CEL-seq","Drop-seq","Microwell-seq"))



seq.col=c("#a5cc26","#A46ACD","#E78D35","#002fa7")


####tool order----------

test = data %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    all.mean.cor = mean(mean_cor, na.rm = TRUE),
    #sd_value = sd(value, na.rm = TRUE),
    .groups = "drop"  
  ) %>%
  arrange(desc(all.mean.cor)) 

data$method<- factor(data$method,
                    levels=test$method)


p=ggplot(data, aes(x = method,y = mean_cor, color =protocol,fill= protocol,group = protocol)) +  
  #geom_line(linetype="dashed",size=1) +
  geom_line() +
  geom_point(size=2,shape=21) + 
  scale_color_manual(values=seq.col)+
  scale_fill_manual(values=seq.col)+
  labs(x = "", y = "correlation coefficient")+p_theme+
  guides(color=guide_legend(title='Protocol'),fill='none')+
  theme(legend.position = "top",
    axis.text.x = element_text(angle = 45,hjust=1),
       
        plot.title = element_text(
          hjust = 0.5,  
          size = 10,    
          face = "bold"  
        ),
        legend.text = element_text(size = 7)) 

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/count.cor/gount.cor.pptx"), width = 5, height = 4)


print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/count.cor/gount.cor.png"),width=5,height=4) 



#########################################################################################


###########二：PAU concordance--------------------------------------------------------------------------

###1.Single cell dataset-----------------------------------------

samples=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
          "arab")
tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","SCAPTURE","Infernape","scUTRquant","scraps","MAAPER")




anno.lists=readRDS("G:/scAPA/summary/temp/fig2/three.anno.lists.rds")
sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
paclists=list()

GSMlists=unique(sample.data$Sample_ID)

for(i in 1: 3){
  #sample=samples[i]
  GSM=GSMlists[i]
  
  
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  ##scPACdataset
  temp.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  
  tool.lists=sample.data[sample.data$Sample_ID==GSM,]$Tool
  
  for(j in 1: length(tool.lists)) {
    tool=tool.lists[j]
    pac=temp.pac[[GSM]][[tool]]
    
    #pac=subsetPACds(pac,totPACtag = ncol(pac@counts)/20,verbose = T)
    
    if(tool=="MAAPER"){
      paclists[[GSM]][[tool]]=pac
      
    } else{
      
      cell.anno=anno.lists[[GSM]]
     
      colnames(pac@counts)=gsub("-.*","",colnames(pac@counts))
      row.names(pac@colData)=gsub("-.*","",row.names(pac@colData))
      pac@colData$barcode=row.names(pac@colData)
      cell.anno$Barcode=gsub("-.*","",cell.anno$Barcode)
      
      loc=match(pac@colData$barcode,cell.anno$Barcode)
      pac@colData$celltype=cell.anno$CellType[loc]
      
      if(anyNA(loc)==TRUE){
        stop(paste0(GSM," (",tool,")",": 细胞类型注释信息不匹配"))
      }
      
      paclists[[GSM]][[tool]]=pac
      
    }
    
    print(paste0(GSM,": ",j,"/",length(tool.lists)))
  
    
  }
}

saveRDS(paclists,file = paste0("G:/scAPA/summary/temp/fig2/bech_pacount/","three.utr.paclists.rds"))


##2.Conservative APA gene ------------------------- 
#identification by at least six methods and present in the reference dataset

info.gene=readRDS("G:/scAPA/summary/temp/fig3/com.gene.rds")


##The difference in pau between proximal pA and distal pA on the gene in single-cell datasets-------------------------------------------
tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","polyApipe","scUTRquant","scraps","MAAPER")
sc_paclists=readRDS(paste0("G:/scAPA/summary/temp/fig2/bech_pacount/","three.utr.paclists.rds"))
#anno.lists=readRDS("G:/scAPA/summary/temp/fig2/three.anno.lists.rds")
info.gene=readRDS("G:/scAPA/summary/temp/fig3/com.gene.rds")


sc_paulists=list()
for(i in 1:3){
  sample=names(sc_paclists)[i]
  tool.lists=names(sc_paclists[[sample]])
  tool.lists=setdiff(tool.lists, "MAAPER")
  
  for(j in 1 :length(tool.lists)){
    message("start")
    
    tool=tool.lists[j]
    pac=sc_paclists[[i]][[tool]]
    
   
    
    if(tool=="MAAPER"){
      ##The MAAPER expression matrix has only one column, which will cause errors later.
      #Therefore, we first add a duplicate column to it.
      
      pac@counts <- cbind(pac@counts, count_copy = pac@counts[, "count"])
      pac@colData<-data.frame(row.names = colnames(pac@counts),
                                group=rep("group1",ncol(pac@counts)))
   
    
    pac=get3UTRAPApd(
      pac,
      minDist = 50,
      maxDist = 10000,
      minRatio = 0.05,
      fixDistal = FALSE,
      addCols = "pd"
    )
    
   srud=subset(pac@anno,pac@anno$pdWhich=="D")
   srud=srud[,c("gene","pdRatio")]
   colnames(srud)=c("gene","rud")
   
   pau <-srud %>%
     mutate(
      
       rup = ifelse(is.nan(rud), 0, 1 - rud),
       rud = ifelse(is.nan(rud), 0, rud)
     )
   row.names(pau)=1:nrow(pau)
   
   pau$pau.diff=pau$rup-pau$rud
   
    }else{
      
      pac=get3UTRAPApd(
        pac,
        minDist = 50,
        maxDist = 10000,
        minRatio = 0.05,
        fixDistal = FALSE,
        addCols = "pd"
      )
      
      
      count <- as.data.frame(as.matrix(pac@counts)) %>%
        tibble::rownames_to_column("paid") %>%
        tidyr::pivot_longer(cols = -paid, names_to = "cell", values_to = "count")
      
      count$celltype <- pac@colData$celltype[match(count$cell, colnames(pac@counts))]
      
   
      count <- count %>%
        group_by(paid, celltype) %>%
        summarise(counts = sum(count)) %>%
        tidyr::pivot_wider(names_from =celltype, values_from = counts) %>%
        as.data.frame()
      
      row.names(count)=count$paid
      count=select(count,-paid)
      ##使@counts和@anno里的行名一致
      count=count[row.names(pac@anno),]
      
      pac@counts=as.matrix(count)
      
      srud=movAPAindex(pac, method="smartRUD", sRUD.oweight=FALSE)
      
      srud=as.data.frame(srud) %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(cols=-gene,names_to = "celltype", values_to = "rud")
      
      pau <-srud %>%
        mutate(
          rup = ifelse(is.nan(rud), 0, 1 - rud),
          rud = ifelse(is.nan(rud), 0, rud)
        )
      pau$pau.diff=pau$rup-pau$rud
      
    }
    
    
   
   ##Keeping the conservative APA gene
   info=info.gene[[sample]]
   pau=subset(pau,pau$gene %in% info)
   
   sc_paulists[[sample]][[tool]]=pau
   
   message(paste0(sample,": ",j,"/",length(tool.lists)))
  }
}

saveRDS(sc_paulists,file="G:/scAPA/summary/temp/fig2/bech_pacount/sc_paulists.rds")
   
##The difference in pau between proximal pA and distal pA on the gene in bulk datasets--------------
bulk_paclists=readRDS("G:/scAPA/summary/temp/fig2/fil.bulk_paclists.rds")


bulk_paulists=list()
for(i in 1:3){
  sample=names(bulk_paclists)[i]
 
    pac=bulk_paclists[[i]]
    
     pac=get3UTRAPApd(
        pac,
        minDist = 50,
        maxDist = 10000,
        minRatio = 0.05,
        fixDistal = FALSE,
        addCols = "pd"
      )
      
  
     
      count <- as.data.frame(as.matrix(pac@counts)) %>%
        tibble::rownames_to_column("paid") %>%
        tidyr::pivot_longer(cols = -paid, names_to = "group", values_to = "count")
      
     
      count$celltype <- pac@colData$group[match(count$group, colnames(pac@counts))]
      
      count$celltype<- gsub("SHAM-FIB", "Sham_Fibroblasts",count$celltype)
      count$celltype<- gsub("SHAM-CD45", "Sham_leukocytes",count$celltype)
      count$celltype<- gsub("SHAM-CD31", "Sham_EC",count$celltype)
      count$celltype<- gsub("MI-CD45", "MI_leukocytes",count$celltype)
      
      
      count <- count %>%
        group_by(paid, celltype) %>%
        summarise(counts = sum(count)) %>%
        tidyr::pivot_wider(names_from =celltype, values_from = counts) %>%
        as.data.frame()
      
      row.names(count)=count$paid
      count=select(count,-paid)
      ##使@counts和@anno里的行名一致
      count=count[row.names(pac@anno),]
      
      pac@counts=as.matrix(count)
      
      srud=movAPAindex(pac, method="smartRUD", sRUD.oweight=FALSE)
      
      srud=as.data.frame(srud) %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(cols=-gene,names_to = "celltype", values_to = "rud")
      
      pau <-srud %>%
        mutate(
          rup = ifelse(is.nan(rud), 0, 1 - rud),
          rud = ifelse(is.nan(rud), 0, rud)
        )
      pau$pau.diff=pau$rup-pau$rud
    
    
    bulk_paulists[[sample]]=pau
    
  }

saveRDS(bulk_paulists,file="G:/scAPA/summary/temp/fig2/bech_pacount/bulk_paulists.rds")



###Using 3seq data as ground truth, calculate the difference in pau across each gene--------------------------------------------


res=data.frame()
for(i in 1:3){
  sample=names(bulk_paulists)[i]
  
  #info=info.gene[[i]]$APA.gene
  if (i == 1) {
    celltypes=c("SC","RS","ES")
  } else if (i == 2) {
    celltypes=c("naive_Tcell","activated_Tcell")
    
  } else if (i == 3) {
    celltypes=c("Sham_Fibroblasts","Sham_leukocytes","Sham_EC","MI_leukocytes")
    
  } else {
    stop("Error：i 必须是 1、2 或 3")
    
  }
  
  
  for(j in 1 :length(tool.lists)){
    message("start")
    bulk_pau=bulk_paulists[[i]]
    
    sc_pau=sc_paulists[[i]]
    tool=tool.lists[j]
    sc_pau=sc_pau[[tool]]
    
    for(c in 1:length(celltypes)){
      
      type=celltypes[c]
      
      bulk_pau2=subset(bulk_pau,bulk_pau$celltype==type)
      
      if(!tool=="MAAPER"){
        sc_pau2 =subset(sc_pau,sc_pau$celltype==type)
      }else{
        sc_pau2 =sc_pau
      }
    
    
    #Using shared genes for comparison
    com=intersect(sc_pau2$gene,bulk_pau2$gene)
    sc_pau2=subset(sc_pau2,sc_pau2$gene %in% com)
    bulk_pau2=subset(bulk_pau2,bulk_pau2$gene %in% com)
    
    
    
    combind<- sc_pau2 %>%
      inner_join(bulk_pau2, by = "gene", suffix = c("_sc", "_bulk"))%>%
      mutate(rud.diff = rud_sc - rud_bulk) %>%
      mutate(rup.diff=rup_sc -rup_bulk) %>%
      mutate(sum.diff=abs(rup.diff)+abs(rud.diff))
    #max(combind$sum.diff)
    
  
    sum.diff=combind[,c("gene","rud.diff","rup.diff","sum.diff")]
    sum.diff$celltype=type
    sum.diff$gene.num=nrow(sum.diff)
    
    # Create cumulative density data
    sum.diff=sum.diff[order(sum.diff$`sum.diff`), ]
    sum.diff$proportion=seq_along(sum.diff$sum.diff) / length(sum.diff$sum.diff)
    
    sum.diff$tool=tool
    sum.diff$sample=sample
    res=rbind(res,sum.diff)
    
    
    message(paste0(sample,"(",tool,")",": ",j,"/",length(tool.lists)))
  }
  }
}

write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/res_sum.diff.csv",row.names = FALSE) 


###the Pearson correlation between the △PAU estimated by different mehtods for each gene and the results from 3’ seq and bulk RNA-seq data--------
bulk_paulists=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/bulk_paulists.rds")
sc_paulists=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/sc_paulists.rds")
tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps")
#info.gene=readRDS("G:/scAPA/summary/temp/fig2/info.gene.rds")

cor.res=data.frame() 
for(i in 1:3){
  sample=names(bulk_paulists)[i]
  
  if (i == 1) {
   celltypes=c("SC","RS","ES")
  } else if (i == 2) {
    celltypes=c("naive_Tcell","activated_Tcell")
    
  } else if (i == 3) {
    celltypes=c("Sham_Fibroblasts","Sham_leukocytes","Sham_EC","MI_leukocytes")
   
  } else {
    stop("Error：i 必须是 1、2 或 3")
    
  }
  
  tool.lists=names(sc_paulists[[i]])
  
  for(j in 1 :length(tool.lists)){
    message("start")
    bulk_pau=bulk_paulists[[i]]
    bulk_pau=as.data.frame(bulk_pau)
    
    sc_pau=sc_paulists[[i]]
    tool=tool.lists[j]
    sc_pau=sc_pau[[tool]]
    sc_pau=as.data.frame(sc_pau)
    
    
   
    for(c in 1:length(celltypes)){
      
      type=celltypes[c]
      
      bulk_pau2=subset(bulk_pau,bulk_pau$celltype==type)
      
      if(!tool=="MAAPER"){
       sc_pau2 =subset(sc_pau,sc_pau$celltype==type)
      }else{
        sc_pau2 =sc_pau
      }
      
      
      if(keep.all.gene==TRUE){
        
        com=intersect(sc_pau2$gene,bulk_pau2$gene)
        sc_pau2=subset(sc_pau2,sc_pau2$gene %in% com)
        bulk_pau2=subset(bulk_pau2,bulk_pau2$gene %in% com) 
        
      }else{  
        
          #Using shared genes
          com=intersect(sc_pau2$gene,bulk_pau2$gene)
          sc_pau2=subset(sc_pau2,sc_pau2$gene %in% com)
          bulk_pau2=subset(bulk_pau2,bulk_pau2$gene %in% com) 
          
          ##The top 60% of genes with smaller △PAU differences
          res.diff <- sc_pau2 %>%
            inner_join(bulk_pau2, by = "gene", suffix = c("_sc", "_bulk")) %>%
            mutate(diff = pau.diff_sc - pau.diff_bulk) %>%
            arrange(abs(diff))  %>% 
            #slice_head(n = 600)
            slice_head(prop= 0.6)
          
          bulk_pau2=subset(bulk_pau2,bulk_pau2$gene %in% res.diff$gene)
          sc_pau2=subset(sc_pau2,sc_pau2$gene %in% res.diff$gene)
        }
        
     
      row.names(bulk_pau2)=bulk_pau2$gene
      row.names(sc_pau2)=sc_pau2$gene
      bulk_pau.diff=select(bulk_pau2,pau.diff)
      sc_pau.diff=select(sc_pau2,pau.diff)
      
      
      cor=psych::corr.test(sc_pau.diff,bulk_pau.diff,  
                           use = "pairwise", 
                           method="pearson", 
                           adjust = "fdr"   
      )
      
      temp=data.frame(sample=sample,
                      method=tool,
                      cor=cor$r[1,1],
                      p.adj=cor$p.adj[1,1],
                      gene.num=cor$n,
                      celltype=type
                      )
    
      cor.res=rbind(cor.res,temp)
      
      
       
    }
    
    message(paste0(sample,": ",j,"/",length(tool.lists)))
    
  }
}

write.csv(cor.res,file="G:/scAPA/summary/temp/fig2/bech_pacount/cor.res.csv",row.names = FALSE)

####Correlation Heatmap-----------------------------------------------------------------
cor.res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cor.res.csv")

#cor.res$method<- factor(cor.res$method,
#                  levels=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))

cor.res$celltype<- factor(cor.res$celltype,
                      levels =c("SC","RS","ES","naive_Tcell","activated_Tcell","Sham_Fibroblasts","Sham_leukocytes","Sham_EC","MI_leukocytes"),
                      labels =c("SC","RS","ES","Naive T cell","Activated T cell","Sham Fibroblasts","Sham leukocytes","Sham EC","MI leukocytes"))

##(1)mouse_sperm##################################################################
data="mouse_sperm"
sub.cor.res=subset(cor.res,cor.res$sample==data)

##tool order(based on gene number)----------------------


test <- sub.cor.res %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean.cor = mean(cor, na.rm = TRUE), 
    mean.gene = mean(gene.num, na.rm = TRUE) 
  ) %>%
  # 添加排名列
  dplyr::mutate(
    cor.rank = rank(-mean.cor, ties.method = "min"),  
    gene.rank = rank(-mean.gene, ties.method = "min"),  
    rank.score = cor.rank *0.3+ gene.rank*0.7,  
    rank = rank(rank.score, ties.method = "min")  
  ) %>%
  arrange(gene.rank)

sub.cor.res$method=factor(sub.cor.res$method,
                          levels = rev(test$method))

col=colorRampPalette(c("#869c6b", "white", "#e58994"), alpha = TRUE)
col2=colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
col3=c("#fee3ce","#eabaa1","#dc917b","#d16d5b","#c44438","#b7282e")



p=ggplot(sub.cor.res, aes(x =celltype , y =method, fill =cor)) + 
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  #scale_fill_gradient2(low = '#3E4F94',
  #                     mid = 'white',
  #                     high = '#b5182b'
  #)+
  scale_fill_gradientn(
    colours = col3)+  #col(10)
  geom_text(aes(label = round(cor,2)), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill="Correlation Coefficient") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank(),         
        axis.title.y = element_blank(),        
        axis.text.y = element_blank(),         
        axis.ticks.y = element_blank()         
  )

##---------------------------------------------------------------------

###每种工具用于计算相关性的基因数目

gene.num=unique(sub.cor.res[,c("method","gene.num")])

sup=ggplot(data=gene.num, 
       aes(x=gene.num, y=method)) +
  geom_bar(stat="identity",fill="#cecccb")+
  p_theme+
  labs(y=NULL,x=bquote("Gene Number"))+
  scale_x_reverse() +  
  scale_y_discrete(position = "right") +  
  guides(fill='none')+
theme(axis.text.x = element_text(angle = 45,hjust=1))
  
  

# 合并图形，对齐 y 轴
p1=plot_grid(
  sup, 
  p, 
  ncol = 2, 
  align = "h",  
  rel_widths = c(0.5, 0.5)  
)


topptx(p1, filename = paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/",data,"_cor.pptx"), width = 6, height = 3)

print(p1)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/",data,"_cor.png"),width=6,height=3)  


##(2)mouse_tcell##################################################################
data="mouse_tcell"
sub.cor.res=subset(cor.res,cor.res$sample==data)


##tool order(based on gene number)----------------------

test <- sub.cor.res %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean.cor = mean(cor, na.rm = TRUE), 
    mean.gene = mean(gene.num, na.rm = TRUE)  
  ) %>%
  dplyr::mutate(
    cor.rank = rank(-mean.cor, ties.method = "min"),  
    gene.rank = rank(-mean.gene, ties.method = "min"),  
    rank.score = cor.rank *0.3+ gene.rank*0.7,  
    rank = rank(rank.score, ties.method = "min")  
  ) %>%
  arrange(gene.rank)

sub.cor.res$method=factor(sub.cor.res$method,
                          levels = rev(test$method))


p=ggplot(sub.cor.res, aes(x =celltype , y =method, fill =cor)) + 
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  
  scale_fill_gradientn(
    colours = col3)+  #col(10)
  geom_text(aes(label = round(cor,2)), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill="Correlation Coefficient") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank(),         
        axis.title.y = element_blank(),        
        axis.text.y = element_blank(),       
        axis.ticks.y = element_blank()         
  )


###gene num

gene.num=unique(sub.cor.res[,c("method","gene.num")])

sup=ggplot(data=gene.num, 
           aes(x=gene.num, y=method)) +
  geom_bar(stat="identity",fill="#cecccb")+
  p_theme+
  labs(y=NULL,x=bquote("Gene Number"))+
  scale_x_reverse() +  
  scale_y_discrete(position = "right") +  
  guides(fill='none')+
  theme(axis.text.x = element_text(angle = 45,hjust=1))



# 合并图形，对齐 y 轴
p2=plot_grid(
  sup, 
  p, 
  ncol = 2, 
  align = "h",  
  rel_widths = c(0.5, 0.5)  
)


topptx(p2, filename = paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/",data,"_cor.pptx"), width = 6, height = 3)

print(p2)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/",data,"_cor.png"),width=6,height=3)  


##(3)mouse_tip##################################################################
data="mouse_tip"
sub.cor.res=subset(cor.res,cor.res$sample==data)

##tool order(based on gene number)----------------------

test <- sub.cor.res %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean.cor = mean(cor, na.rm = TRUE), 
    mean.gene = mean(gene.num, na.rm = TRUE)  
  ) %>%
  dplyr::mutate(
    cor.rank = rank(-mean.cor, ties.method = "min"), 
    gene.rank = rank(-mean.gene, ties.method = "min"),  
    rank.score = cor.rank *0.3+ gene.rank*0.7, 
    rank = rank(rank.score, ties.method = "min")  
  ) %>%
  arrange(gene.rank)

sub.cor.res$method=factor(sub.cor.res$method,
                          levels = rev(test$method))



p=ggplot(sub.cor.res, aes(x =celltype , y =method, fill =cor)) + 
  geom_tile(
    color = "#FAFAFA", 
    lwd = 0.8, 
    linetype = 1) + 
  scale_fill_gradientn(
    colours = col3)+  #col(10)
  geom_text(aes(label = round(cor,2)), 
            size = 4, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),fill="Correlation Coefficient") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        legend.title = element_text(size = 7), 
        panel.grid = element_blank(),        
        axis.title.y = element_blank(),       
        axis.text.y = element_blank(),        
        axis.ticks.y = element_blank()         
  )



###gene number

gene.num=unique(sub.cor.res[,c("method","gene.num")])

sup=ggplot(data=gene.num, 
           aes(x=gene.num, y=method)) +
  geom_bar(stat="identity",fill="#cecccb")+
  p_theme+
  labs(y=NULL,x=bquote("Gene Number"))+
  scale_x_reverse() +  
  scale_y_discrete(position = "right") +  
  guides(fill='none')+
  theme(axis.text.x = element_text(angle = 45,hjust=1))



# 合并图形，对齐 y 轴
p3=plot_grid(
  sup, 
  p, 
  ncol = 2, 
  align = "h", 
  rel_widths = c(0.5, 0.5)  
)


topptx(p3, filename = paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/",data,"_cor.pptx"), width = 6, height = 3)

print(p3)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/",data,"_cor.png"),width=6,height=3)  



####Merge figures#####################################################################
p=plot_grid(p2,p1,p3,nrow = 3)

topptx(p, filename = paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/cor.pptx"), width = 6, height = 9)

print(p)
graph2png(file=paste0("E:/scAPA/summary/Figure/fig2/bech_pacount/cor.png"),width=6,height=9) 





#############三：Cell Clustering--------------------------------------------------------

###1.Gene expression matrix-------------------------------------------------------
GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
           "GSM2333586","GSM2333587", "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
           "GSM3490689","GSM3490690","GSM3490691")


sub.GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
               "3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490690","GSM3490691")


sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
#anno.lists=readRDS("G:/scAPA/summary/temp/fig2/three.anno.lists.rds")
anno.lists=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/anno.lists.rds")
gcounts=readRDS("G:/scAPA/summary/gene/gcounts.rds")



res=data.frame()
obj.lists=list()

for(i in 1:length(sub.GSMlists)){
  GSM=sub.GSMlists[i]
  counts=gcounts[[GSM]]
  
  cell.anno=anno.lists[[GSM]]
  cell.anno$Barcode=gsub("-.*","",cell.anno$Barcode)
  colnames(counts)=gsub("-.*","",colnames(counts))
  
  info=intersect(colnames(counts),cell.anno$Barcode)
  counts=counts[,info]
  cell.anno=subset(cell.anno,cell.anno$Barcode %in% info)
  
  loc=match(colnames(counts),cell.anno$Barcode)
  if(anyNA(loc)==TRUE){
    stop(paste0(GSM,": 细胞类型注释信息不匹配"))
  }
  
  
  anno=data.frame(row.names = colnames(counts),
                  barcode=colnames(counts),
                  celltype=cell.anno$CellType[loc])
  
  #keep genes expressed in at least 0.5% of cells
  obj<- CreateSeuratObject(counts =counts,
                           meta.data = anno,
                           min.cells = 2,
                           min.features =ncol(counts)/200, 
                           project = "gcounts" )
  
  obj<- NormalizeData(obj)
  ##Number of high-variability features
  n_features <- length(rownames(obj))
  if( n_features >= 3000){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  } else if(n_features >= 2000){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
  } else if(n_features >= 1000){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 800, verbose = FALSE)
  } else if(n_features >= 500){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 500, verbose = FALSE)
  } else {
   
    VariableFeatures(obj) <- rownames(obj)
    message("Gene number <500: Using all genes as variable features.")
  }
  
  obj<- ScaleData(obj, verbose = FALSE)
  
  obj <- RunPCA(obj , npcs = 10, verbose = FALSE)
  
  
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:10)
  obj <- FindNeighbors(obj, dims = 1:10)
  
  n_types <- length(unique(obj@meta.data$celltype))  
  # Dynamically set the resolution
  if(n_types <= 4){
    resolution <- 0.1
  }else if(n_types < 8){
    resolution <- 0.2
  }else if(n_types <= 10){
    resolution <- 0.3
  } else if(n_types <= 15){
    resolution <- 0.5
  } else if(n_types <= 20){
    resolution <- 0.7
  } else {
    resolution <- 1.0
  }
  obj <- FindClusters(obj, resolution =resolution)
  
  obj.lists[[GSM]]=obj
  
  ##Clustering results and real cell labels
  true_labels=obj@meta.data$celltype
  cluster_labels=obj@meta.data$seurat_clusters
  
  ### External metrics----------------------------------------------------------
  if(!is.null(true_labels)){
    
    ari <- aricode::ARI(true_labels, cluster_labels)
    #ari <- mclust::adjustedRandIndex(true_labels, cluster_labels)
    

    nmi <- aricode::NMI(true_labels, cluster_labels)
    print(paste("ARI:", ari, "NMI:", nmi))
  }
  
  
  ###Internal metrics------------------------------------------------------------------
  cluster_labels<- as.numeric(as.character(cluster_labels)) + 1
  
  if (length(unique(cluster_labels)) == 1) {
    sc <- 0
    db <- 0
    message(paste0(tool,": 只有一种聚类标签，无法计算SC，DBI"))
  } else{
    
    data_matrix <- GetAssayData(obj, slot = "scale.data")
    dist_matrix <- dist(t(data_matrix))
    
   
    sil <- silhouette(as.numeric(cluster_labels), dist_matrix)
    sc <- mean(sil[, "sil_width"])
    print(paste("Mean Silhouette Width:", sc))
    
    
    # Davies-Bouldin
    db <- clusterSim::index.DB(t(data_matrix), as.numeric(cluster_labels))$DB
    print(paste("Davies-Bouldin Index:", db))
    
    
  }
  
  
    temp=data.frame(sample=GSM,
                    tool="CellRanger/StarSolo",
                    method=c("ARI","NMI","SC","DBI"),
                    value=paste0(c(ari,nmi,sc,db)),
                    gene.number=n_features,
                    var.gene.number=length(VariableFeatures(obj))
    )
    
    res=rbind(res,temp) 
  
  
  
  
  message(paste0(GSM,": ",i,"/",length(sub.GSMlists)))
}

res$value <- as.numeric(res$value)
res$type="gcounts"

sum(is.na(res$value))

write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster.gcounts.res.csv",row.names = FALSE)

saveRDS(obj.lists,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster.gcounts.obj.rds")

res=res[res$method=="ARI",] %>%
  arrange(desc(value))


------------------------------------------------------------------
    
###2.APA data-----------------------------------------------------

sub.GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
               "3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490690","GSM3490691")

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
#anno.lists=readRDS("G:/scAPA/summary/temp/fig2/three.anno.lists.rds")
anno.lists=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/anno.lists.rds")


#info.gene=readRDS("G:/scAPA/summary/temp/fig2/info.gene.rds")

###Select the PAC data type for clustering------------------

keep.all.gene=TRUE
data="rud"  #srud/rud /count

#########################################


res=data.frame()
obj.lists=list()

#for(i in 10:length(sub.GSMlists)){
for(i in 9:17){
  GSM <-  sub.GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  all.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  
  tool.lists=names(all.pac[[GSM]])
  tool.lists=setdiff(tool.lists, "MAAPER")
  
  
  for(j in 1 :length(tool.lists)){
    #message("start")
    
    tool=tool.lists[j]
    pac=all.pac[[GSM]][[tool]]
    
    
    ##Should the conservative APA gene be used?
    if(keep.all.gene==FALSE){
      info=info.gene[[i]]$APA.gene
      pac@anno=subset(pac@anno,pac@anno$gene %in% info)
      pac@counts=pac@counts[row.names(pac@anno),]
    }
    
    
    ##Add cell type annotation information
    
    cell.anno=anno.lists[[GSM]]
    
    colnames(pac@counts)=gsub("-.*","",colnames(pac@counts))
    row.names(pac@colData)=gsub("-.*","",row.names(pac@colData))
    pac@colData$barcode=row.names(pac@colData)
    cell.anno$Barcode=gsub("-.*","",cell.anno$Barcode)
    
    if((tool=="SCAPE" & GSM=="GSM3490690") | (tool=="SCAPE" & GSM=="GSM3490691")){
      colnames(pac@counts)=sub("^.*\\.", "",colnames(pac@counts)) 
      row.names(pac@colData)=sub("^.*\\.", "",row.names(pac@colData))
      pac@colData$barcode=row.names(pac@colData)
      
    }
    
    
    loc=match(pac@colData$barcode,cell.anno$Barcode)
    if(anyNA(loc)==TRUE){
      stop(paste0(GSM,": 细胞类型注释信息不匹配"))
    }
    
    
    pac@colData$celltype=cell.anno$CellType[loc]
    
    #####SmartRU index
    if(data=="srud"){
      pac=get3UTRAPApd(
        pac,
        minDist = 50,
        maxDist = 10000,
        minRatio = 0.05,
        fixDistal = FALSE,
        addCols = "pd"
      )
      
      srud=movAPAindex(pac, method="smartRUD", sRUD.oweight=FALSE)
      #table(is.nan(srud))
      srud[is.nan(srud)]=0
      
      obj<- CreateSeuratObject(counts =srud,
                               meta.data = pac@colData,
                               min.cells = 0, 
                               project = "RUD" )
      
      
    } else if(data=="rud"){
      ds=get3UTRAPAds(pac, sortPA=TRUE, choose2PA=NULL)
      #summary(ds)
      
      ##RUD index
      rud=movAPAindex(ds, method='rud', choose2PA=NULL)
      #table(is.nan(srud))
      rud[is.nan(rud)]=0
      
      obj<- CreateSeuratObject(counts =rud,
                               meta.data = pac@colData,
                               min.cells = 0, 
                               project = "RUD" )
      
      
    } else if(data=="rup"){
      ds=get3UTRAPAds(pac, sortPA=TRUE, choose2PA=NULL)
      #summary(ds)
      
      ##RUP index
      rup=movAPAindex(ds, method='rup', choose2PA=NULL)
      #table(is.nan(srud))
      rup[is.nan(rup)]=0
      
      obj<- CreateSeuratObject(counts =rup,
                               meta.data = pac@colData,
                               min.cells = 0, 
                               project = "RUP" )
      
      
    }  else if(data=="psi"){
      ds=get3UTRAPAds(pac, sortPA=TRUE, choose2PA=NULL)
      #summary(ds)
      
      ##PSI
      psi=movAPAindex(ds, method='psi', choose2PA=NULL)
      #table(is.nan(srud))
      psi[is.nan(psi)]=0
      
      obj<- CreateSeuratObject(counts =psi,
                               meta.data = pac@colData,
                               min.cells = 0, 
                               project = "PSI" )
      
      
    } else if(data=="count"){
      
      obj<- CreateSeuratObject(counts =pac@counts,
                               meta.data = pac@colData,
                               min.cells = 0, 
                               project = "PA count" )
      
    }else {
      stop("Error：data must be  pA ratio or pA count")
      
    }
    
    
    obj<- NormalizeData(obj)
   
    n_features <- length(rownames(obj))
    if( n_features >= 3000){
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    } else if(n_features >= 2000){
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
    } else if(n_features >= 1000){
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 800, verbose = FALSE)
    } else if(n_features >= 500){
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 500, verbose = FALSE)
    } else {
      
      VariableFeatures(obj) <- rownames(obj)
      message("Gene number <500: Using all genes as variable features.")
    }
    
    obj<- ScaleData(obj, verbose = FALSE)
    
    obj <- RunPCA(obj , npcs = 10, verbose = FALSE)
    
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:10)
    obj <- FindNeighbors(obj, dims = 1:10)
    
    n_types <- length(unique(obj@meta.data$celltype))  
    # Dynamically set the resolution.
    if(n_types <= 4){
      resolution <- 0.1
    }else if(n_types < 8){
      resolution <- 0.2
    }else if(n_types <= 10){
      resolution <- 0.3
    } else if(n_types <= 15){
      resolution <- 0.5
    } else if(n_types <= 20){
      resolution <- 0.7
    } else {
      resolution <- 1.0
    }
    
    obj <- FindClusters(obj, resolution =resolution)
    
    obj.lists[[GSM]][[tool]]=obj
    
    ##Clustering results and true cell labels
    true_labels=obj@meta.data$celltype
    cluster_labels=obj@meta.data$seurat_clusters
    
    ### External metrics----------------------------------------------------------
    if(!is.null(true_labels)){
      
     
      ari <- aricode::ARI(true_labels, cluster_labels)
      #ari <- mclust::adjustedRandIndex(true_labels, cluster_labels)
      
     
      nmi <- aricode::NMI(true_labels, cluster_labels)
      print(paste("ARI:", ari, "NMI:", nmi))
    }
    
    
    ###Internal metrics------------------------------------------------------------------
   
    cluster_labels<- as.numeric(as.character(cluster_labels)) + 1
    
    if (length(unique(cluster_labels)) == 1) {
      sc <- 0
      db <- 0
      message(paste0(tool,": 只有一种聚类标签，无法计算SC，DBI"))
    } else{
      
      data_matrix <- GetAssayData(obj, slot = "scale.data")
      dist_matrix <- dist(t(data_matrix))
      
      
      sil <- silhouette(as.numeric(cluster_labels), dist_matrix)
      sc <- mean(sil[, "sil_width"])
      print(paste("Mean Silhouette Width:", sc))
      
      
      # Davies-Bouldin
      db <- clusterSim::index.DB(t(data_matrix), as.numeric(cluster_labels))$DB
      print(paste("Davies-Bouldin Index:", db))
      
    }
    
    if(data=="count"){
      
      temp=data.frame(sample=GSM,
                      tool=tool,
                      method=c("ARI","NMI","SC","DBI"),
                      value=paste0(c(ari,nmi,sc,db)),
                      pA.number=n_features,
                      var.pA.number=length(VariableFeatures(obj))
      )
      
      res=rbind(res,temp)
    } else{
      temp=data.frame(sample=GSM,
                      tool=tool,
                      method=c("ARI","NMI","SC","DBI"),
                      value=paste0(c(ari,nmi,sc,db)),
                      gene.number=n_features,
                      var.gene.number=length(VariableFeatures(obj))
      )
      
      res=rbind(res,temp) 
    }
    
    
    
    message(paste0(GSM,": ",j,"/",length(tool.lists)))
  }
}


res$value <- as.numeric(res$value)
res$type=data

sum(is.na(res$value))

write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.rud.res(human+arab).csv",row.names = FALSE)
#write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster.res(info.gene).csv",row.names = FALSE)
saveRDS(obj.lists,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.rud.obj(human+arab).rds")


###(1)Summarized clustering metric scores----------------------------
##Merge data of the same type
res1=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.rud.res(mouse).csv")
res2=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.rud.res(human+arab).csv")

res=rbind(res1,res2)
write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.rud.res.csv",row.names = FALSE)

##Merge the results of the three types of data（gcounts, pa count,rud)
g.res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.gcounts.res.csv")
colnames(g.res)[c(5,6)]=c("feature.num","var.feature.num")

res1=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.count.res.csv")
colnames(res1)[c(5,6)]=c("feature.num","var.feature.num")

res2=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.rud.res.csv")
colnames(res2)[c(5,6)]=c("feature.num","var.feature.num")

res=rbind(g.res,res1,res2)

write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.res.csv",row.names = FALSE)

###(2)plot------------------------
index.type="ARI"

data=subset(res,res$method==index.type)
 
  
  ##tool order
 order= res[(res$method==index.type & !res$type=="rud"),] %>%
   group_by(tool,type)%>%
   summarise(mean.value=mean(value, na.rm = TRUE)
             )%>%
   arrange(desc(mean.value))
  
 data$tool=factor(data$tool,
                   levels = c("CellRanger/StarSolo",setdiff(order$tool, "CellRanger/StarSolo")),
                   labels = c("Gene Counts",setdiff(order$tool, "CellRanger/StarSolo"))
)

data$type=factor(data$type,
                 levels = c("gcounts","count", "rud"),
                 labels = c("Gene Counts","pA Counts", "RUD"))
cols=c("Gene Counts"="#fff9c4",
      "pA Counts"="#6acfc7",
      "RUD"="#e0f7e6")

cols=c("Gene Counts"="#E6CDC6",
       "pA Counts"="#B3B69C",
       "RUD"="#a0afb6")


###(2.1)barplot--------------------------
y.lab="ARI"

p=ggplot(data, aes(x = tool, y = round(value, 2), fill = type)) +
  geom_bar(
    stat = "summary", 
    fun = "mean",     
    position = position_dodge(width = 0.8),  
    width = 0.7,       
    color = "black",   
    alpha = 0.8        
  ) +
  
  stat_summary(
    geom = "errorbar",
    fun.data = mean_se,  
    position = position_dodge(width = 0.8),
    width = 0.2,         
    color = "black"
  ) +
  scale_fill_manual(values = cols) + 
  labs(
    x = "", 
    y = y.lab, 
    title = "" ,
    fill = ""  
  ) +
  p_theme +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "top" 
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/",index.type,".barplot2.pptx"), width = 6, height = 5)

print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/",index.type,".barplot2.png"),width=6,height=5) 


###(2.2)heatmap---------------------------------------------
res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.res.csv")
sample.data<- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
 loc=match(res$sample,sample.data$Sample_ID)
res$species=sample.data$Species[loc]
res$protocol=sample.data$Protocol[loc] 
 
index.type="ARI"
data.type="rud"



data=subset(res,res$method==index.type & res$type==data.type)

data$sample=factor(data$sample,
       levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
                "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","GSM2333586","GSM2333587","GSM2906396","GSM2906399",
                "3k","4k","5k","6k","GSM4712885","GSM4712907",
                "GSM3490689","GSM3490690","GSM3490691"),
       labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell","Mouse_GFP1","Mouse_GFP2",
                  "Mouse_Intestine1","Mouse_Intestine2","Mouse_Intestine3","Mouse_Intestine4","Mouse_Intestine5","Mouse_Intestine6",
                  "Mouse_Hypothalamus1", "Mouse_Hypothalamus2","Mouse_Bone1","Mouse_Bone2",
                  "Human_Pbmc1","Human_Pbmc2","Human_Pbmc3","Human_Pbmc4","Human_Covid_Pbmc1","Human_Covid_Pbmc2",
                  "Arabidopsis1","Arabidopsis2","Arabidopsis3"))



##tool order based on average ARI score
order= res[(res$method==index.type & res$type==data.type),] %>%
  group_by(tool,type)%>%
  summarise(mean.value=mean(value, na.rm = TRUE)
  )%>%
  arrange(desc(mean.value))

data$tool=factor(data$tool,
                 levels = c("CellRanger/StarSolo",setdiff(order$tool, "CellRanger/StarSolo")),
                 labels = c("Gene Counts",setdiff(order$tool, "CellRanger/StarSolo"))
)

data$type=factor(data$type,
                 levels = c("gcounts","count", "rud"),
                 labels = c("Gene Counts","pA Counts", "RUD"))

data=subset(data,!data$tool=="Gene Counts")

#col <- colorRampPalette(c("#108b96","white","#f3993a"), alpha = TRUE)
#col <- colorRampPalette(c("#2b2e77","white","#f29a76"), alpha = TRUE)

if(data.type=="count"){
  size.lab="pA Number"
}else if(data.type=="rud"){
  size.lab="APA Gene Number"
}

if(index.type=="ARI"){
  legend.lab="ARI"
  col <- colorRampPalette(c("#108b96","white","#f3993a"), alpha = TRUE)
  
}else if(index.type=="SC"){
  legend.lab="Silhouette Coefficient"
  col <- colorRampPalette(c("#2b2e77","white","#f29a76"), alpha = TRUE)
}

p=ggplot(data, aes(x =tool, y =sample,fill=value)) + #fill =-1*log10(p.adj)
  geom_tile(aes(fill = rectheat),color="black",fill="white") + 
  geom_point(
    aes(size = feature.num,         
        fill =value          
        #,shape = is_significant  
        ), 
    shape=21,
   color = "#909fac",     
    stroke = 0.8         
  ) +

  scale_size_continuous(
  range = c(2, 7) ,    
    guide = guide_legend(         
    override.aes = list(shape = 21)))+
scale_fill_gradientn(
  colors =col(8))+  #col2(10)
  geom_text(aes(label = round(value, 2)), 
            size = 3, 
            color = "black") + 
  theme_minimal()+
  theme_bw()+
  labs(x=element_blank(),y=element_blank(),
       title = "",
       fill=legend.lab,
       size=size.lab) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        axis.text.y = element_text(angle = 45, hjust = 1,size = 7),
        panel.background = element_rect(fill =  "gray100", color = 'white'),
        panel.grid.major = element_line(color = 'white', linetype = 'dotted'),
        panel.grid.minor = element_line(color = 'green', size = 2),
        legend.title = element_text(size = 7),  
        panel.grid = element_blank(),        
        axis.title.y = element_blank(),        
        #axis.text.y = element_blank(),         
        axis.ticks.y = element_blank()        
  ) 


topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/",index.type,"_",data.type,".pptx"), width = 6, height = 5)

print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/",index.type,"_",data.type,".png"),width=6,height=5) 

###(2.2.1)Add legend---------------------------------------------------------------------------
temp=data
temp$sub.sample=paste0(temp$species,"_",temp$protocol)
temp$sub.sample= factor(temp$sub.sample,
                              levels =c("Mouse_10x","Mouse_CEL-seq","Mouse_Drop-seq","Mouse_Microwell-seq","Human_10x","Arabidopsis_10x" ),
                              labels = c("Mouse (10X)","Mouse (CEL-seq)","Mouse (Drop-seq)","Mouse (Microwell-seq)","Human (10X)","Arabidopsis (10X)"))

temp$anno <- "Sample_ID"

p1.sample <-ggplot(  temp, aes(y=sample, x=anno, fill=sample)) + geom_tile() + 
  scale_y_discrete(position="left") +
  scale_fill_manual(breaks = temp$sample,
                    values=sample.col,
                    name="Sample_ID")+
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        #axis.text.x = element_text(angle = 45,hjust=1),
        axis.ticks.x = element_blank(),
        axis.text.y=element_blank())+
  xlab(NULL) + ylab(NULL) +guides(fill="none")




temp$anno2 <- "sub.sample"


p1.special <- ggplot(temp, aes(y=sample, x=anno2, fill=sub.sample)) + geom_tile() + 
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

#p$p %>% insert_left(p$p1,width=.1) %>% insert_left(p$p2,width=.1)

sup.p<-(p + theme(axis.text.y = element_blank())) %>% 
  aplot::insert_left((legend.data$one.sample), width=.1)  %>% 
  aplot::insert_left((legend.data$special), width=.1)

topptx(sup.p, filename = paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/",index.type,"_",data.type,"2.pptx"), width = 6, height = 5)

print(sup.p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/",index.type,"_",data.type,"2.png"),width=6,height=5) 



######(3)umap plot---------------------------------------------------------------------
#聚类指标得分
res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/cluster.res.csv")

data="rud"
file.path="G:/scAPA/summary/temp/fig2/bech_pacount/cluster"

obj.lists=readRDS(paste0(file.path,"/cluster.",data,".obj(1-3).rds"))

##tool order based on average ARI score
order= res[(res$method=="ARI" & !res$type=="rud"),] %>%
  group_by(tool,type)%>%
  summarise(mean.value=mean(value, na.rm = TRUE)
  )%>%
  arrange(desc(mean.value))


tool.lists2=c(setdiff(order$tool, "CellRanger/StarSolo"))

col.list=list()
col.list[["mouse_sperm"]]=c('#002fa7',"#387a88","#f9cb45")
col.list[["mouse_tcell"]]=c("#387a88","#f9cb45")
col.list[["mouse_tip"]]=c('#002fa7',"#387a88","#f9cb45","#b3a9eb")
#cluster.col=c('#a5cc26',"#0f86a9","#fc8452","#ff7bac")

plot=list()
for(i in 1:length(obj.lists)){
  GSM=names(obj.lists)[i]
  
 for(j in 1: length(tool.lists2)){
   tool=tool.lists2[j]
   obj=obj.lists[[i]][[j]]
   
  
   p=plotUMAP(obj,group="celltype",xend_add=2,yend_add=4,annotate = FALSE,allcolour =col.list[[i]], label_size=2,point.size=0.1,label = FALSE)+
     labs(title =tool) +  
     theme(plot.title = element_text(size = 15,
                                     face="bold"),
                   plot.background = element_rect(
             colour = "black",  
             fill = NA,        
            size = 0.8 ))+         
     theme(legend.position = "none") 
   
  
                                    
   sup=DimPlot(obj, reduction = "umap", group.by ="celltype", pt.size = 0.1, label = FALSE, label.size = 5, repel = TRUE, shuffle = TRUE)+
     ggplot2::scale_color_manual(values =col.list[[i]])+
     theme(legend.text = element_text(size = 7))
  
   legends <- cowplot::get_plot_component(sup, "guide-box", return_all = TRUE)
   
  
   combined_legend <- plot_grid(plotlist = legends, ncol = 1)
   
   
   
   plot[[GSM]][[tool]]=p
   plot[[GSM]][["legend"]]=combined_legend
   
   message(paste0(GSM,": ",j,"/",length(tool.lists2)))
 
   }
  
}

saveRDS(plot,file=paste0("G:/scAPA/summary/temp/fig2/bech_pacount/cluster/plot.umnp_",data,"(1-3).rds"))



###Arrange by column
cowplot::plot_grid(plotlist = plot[[3]][tool.lists2],  
                   ncol = 2,
                   align = "hv",  
                   axis ="tblr" ) 



sample="mouse_tip"


graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/UMAP/",sample,"_",data,"(col).png"),width=6,height=10)



####四：The correlation of expression levels of consensus pAs-------------------------------------------

###(1.1)Pearson correlation coefficients between pA expression profiles of consensus pAs quantified by different methods-------------------------
#Consensus pAs are those commonly identified by all the ten methods.

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)



com.pa=read.csv("G:/scAPA/summary/temp/fig1.2/com.pa.stat.csv")
com.pa<- com.pa %>%
  group_by(sample_id) %>%
  filter(n_tools_detected == max(n_tools_detected))


##datasets with > 200 consensus pAs were used
sub.GSMlists=c("GSM3629847+GSM3629848","GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP",
               "MI_DAY7_GFP","SHAM_DAY7_GFP",
               "3k","4k","5k","6k","GSM4712885","GSM4712907")
table(com.pa$sample_id)


com.counts=list()
for(i in 1:length( sub.GSMlists)) {
  #for(id in 15:15) {
  GSM <-  sub.GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  all.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  paclists=all.pac[[GSM]]
  tool.lists=names(paclists)
  
  for(j in 1: length(tool.lists)){
    tool=tool.lists[j]
    
    if(tool=="MAAPER"){
      paclists[[j]]@colData$group=tool
    }else{
      paclists[[j]]@colData$barcode=names(all.pac[[GSM]])[j]
    }
    
    
  }
  
  new.pac=mergePACds(paclists, d = 24, by = "coord")
  #head(new.pac@counts)
  #head(new.pac@anno)
  
  new.pac@anno$paid=paste0(new.pac@anno$chr,":",new.pac@anno$strand,":",new.pac@anno$coord)
  new.pac@counts=new.pac@counts[row.names(new.pac@anno),]
  row.names(new.pac@anno)=new.pac@anno$paid
  row.names(new.pac@counts)=row.names(new.pac@anno)
  
  info=com.pa[com.pa$sample_id==GSM,]$paid
  print(paste0(GSM," : ",length(info)," same pA "))
  
  new.pac=subsetPACds(new.pac,PAs =info,verbose =TRUE)
  
  
  
  res <- data.frame(row.names =row.names(new.pac@anno))
  for(j in 1: length(tool.lists)){
    tool=tool.lists[j]
    temp.pac=subsetPACds(new.pac,group = "group",conds = tool)
    temp=data.frame(row.names = row.names(temp.pac@anno),
                    counts=Matrix::rowSums(temp.pac@counts))
    colnames(temp)=tool
    
    res <- cbind(res, temp)  
    
    
  }
  
  com.counts[["counts"]][[GSM]]=res
  
  #Calculate the correlation of expression levels
  cor.test=psych::corr.test(res,method = "pearson",adjust="BH")	
  
  cor=as.data.frame(cor.test[["r"]])
  padj=as.data.frame(cor.test[["p.adj"]])
  
  com.counts[["cor"]][[GSM]]=cor
  com.counts[["padj"]][[GSM]]=padj
  
  message(paste0(GSM,": ",i,"/",length(sub.GSMlists)))
}

saveRDS(com.counts,file="G:/scAPA/summary/temp/fig2/bech_pacount/com.counts.rds")


##(2.2)The expression correlation of consensus pAs between pairwise tools-------------------

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)


data.lists=readRDS("G:/scAPA/summary/temp/fig1.2/new.anno.paclists.rds")


GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
  "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","GSM2333586","GSM2333587","GSM2906396","GSM2906399",
  "3k","4k","5k","6k","GSM4712885","GSM4712907",
  "GSM3490689","GSM3490690","GSM3490691")




pair.com.counts=list()
res=data.frame()
cor.res=data.frame()
for(i in 1:length( GSMlists)) {
  #for(id in 15:15) {
  GSM <-  GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  
  all.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  paclists=all.pac[[GSM]]
  tool.lists=names(paclists)
  
  
  for(j in 1: length(tool.lists)){
    tool=tool.lists[j]
    
    if(tool=="MAAPER"){
      paclists[[j]]@colData$group=tool
    }else{
      paclists[[j]]@colData$barcode=names(all.pac[[GSM]])[j]
    }}
    
   if(GSM=="GSM3629847+GSM3629848"){
     new.pac=readRDS("G:/scAPA/summary/PACds/mouse_esc_merge.pacds.rds")
     
   }else{
     new.pac=mergePACds(paclists, d = 24, by = "coord")
     #head(new.pac@counts)
     #head(new.pac@anno)
     
     new.pac@anno$paid=paste0(new.pac@anno$chr,":",new.pac@anno$strand,":",new.pac@anno$coord)
     new.pac@counts=new.pac@counts[row.names(new.pac@anno),]
     row.names(new.pac@anno)=new.pac@anno$paid
     row.names(new.pac@counts)=row.names(new.pac@anno)
     
   }
 
  
  for(j in 1: length(tool.lists)){
    tool=tool.lists[j]
    
    
  count.data<- data.frame()
 #1.Count the commonly identified pA sites between each pair of tools
  for(k in tool.lists){
   
    data=data.lists[[GSM]]
   
    pa.num <- sum(data@counts[,tool] != 0 & !is.na(data@counts[,tool]))
    temp.pe.num <-sum(data@counts[,k] != 0 & !is.na(data@counts[,k]))
    
   
    info.pa <- sum(
      data@counts[,tool] != 0 & !is.na(data@counts[,tool]) & 
        data@counts[,k] != 0 & !is.na(data@counts[,k])
    )
    
  
   union.pa <- sum(
     (data@counts[,tool] != 0 & !is.na(data@counts[,tool])) | 
      ( data@counts[,k] != 0 & !is.na(data@counts[,k]))
   )
    
    
  
    temp <- data.frame(row.names = NULL,
                       species=species,
                       protocol=protocol,
                       sample=sample,
                       sample_id=GSM,
                       method = tool, 
                       compare.method=k,
                       pa.num=pa.num,
                       compare.pa.num= temp.pe.num,
                       type=c("common","union"),
                       number=c(info.pa,union.pa),
                       
                       jaccard.index=info.pa/(union.pa)
    )
    
   res=rbind(res,temp)

  ##2.Calculate the total expression of pAs commonly identified by each pair of tools
    keep <- data@counts[, tool] != 0 & !is.na(data@counts[, tool]) & 
      data@counts[, k] != 0 & !is.na(data@counts[, k])
    
     keep=row.names(data@counts)[keep]

    print(paste0(GSM," (",tool,"_",k,") ", ": ",length( keep)," same pA "))
    
    sub.new.pac=subsetPACds(new.pac,PAs =keep,verbose =TRUE)
    
 
      temp.pac=subsetPACds( sub.new.pac,group = "group",conds = c(tool))
      temp=data.frame(row.names = row.names(temp.pac@anno),
                      counts=Matrix::rowSums(temp.pac@counts))
      colnames(temp)=tool
      #colnames(temp)="method"
      
      temp.pac2=subsetPACds( sub.new.pac,group = "group",conds = c(k))
      temp2=data.frame(row.names = row.names(temp.pac2@anno),
                      counts=Matrix::rowSums(temp.pac2@counts))
      #colnames(temp2)=k
      colnames(temp2)="compare.method"
      
      count.temp <- cbind(temp, temp2) 
      count.temp$group=paste0(tool,"_",k)
      
      
      count.data=rbind(count.data,count.temp)
      
      
      ##3.Calculate the correlation of expression levels
      #cor=cor(res,method = "pearson") 
      #cor.test=psych::corr.test(res,method = "pearson",adjust="BH")	
      cor=psych::corr.test(count.temp[,1],  count.temp[,2],
                           use = "pairwise", 
                           method="pearson", 
                           adjust = "fdr"   
      )
      
      cor.temp <- data.frame(row.names = NULL,
                         species=species,
                         protocol=protocol,
                         sample=sample,
                         sample_id=GSM,
                         method = tool, 
                         compare.method=k,
                        cor=cor$r,
                        padj=cor$p.adj
      )
      
     cor.res=rbind(cor.res,cor.temp)
      
  }
  
  pair.com.counts[[GSM]][[tool]]= count.data
  message(paste0(GSM,"(",tool,")",": ",i,"/",length(GSMlists)))
  
  }
    
}
saveRDS(pair.com.counts,file="G:/scAPA/summary/temp/fig2/bech_pacount/pair.com.counts.rds")
write.csv(res,file="G:/scAPA/summary/temp/fig2/bech_pacount/pair.jaccard.index.csv",row.names = FALSE)
write.csv(cor.res,file="G:/scAPA/summary/temp/fig2/bech_pacount/pair.count.cor.csv",row.names = FALSE)



##3.plot-----------------------
mycol <- colorRampPalette(c("#cdd5ae", "white", "#e58994"), alpha = TRUE) #"#cdd5ae"
mycol2 <- colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)


##(3.1)Consensus pAs identified by all the ten methods----------------------
com.counts=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/com.counts.rds")
sub.GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
               "MI_DAY7_GFP","SHAM_DAY7_GFP",
               "3k","4k","5k","6k","GSM4712885","GSM4712907")

res=data.frame()
for(i in 1:length( sub.GSMlists)) {
  #for(id in 15:15) {
  GSM <-  sub.GSMlists[i]
  
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
  
  data=com.counts[["cor"]][[GSM]]
  data$method=row.names(data)
  #将宽数据转换为长数据
 data=data%>%
    pivot_longer(
      cols = -method,  
      names_to = "compare.method",  
      values_to = "cor"  
    ) 
 data=as.data.frame(data)

data$sample_id=GSM
data$sampe=sample
data$protocol=protocol
data$species=species

res=rbind(res,data)

}

input.data= res%>%
  group_by(method, compare.method) %>%
  summarise(
    mean_cor = mean(cor),
    num_samples = n(),
    .groups = 'drop'
  )

#转换成宽数据格式，即相关性矩阵
input.data=tidyr::pivot_wider(input.data[,c(1:3)], 
                       names_from = compare.method, 
                       values_from = mean_cor)
  
  
input.data=as.data.frame(input.data)  
row.names(input.data)=input.data$method
input.data=input.data[,-1]
  
p=ggcorrplot::ggcorrplot(
  input.data,
  method = "circle",
  type = "upper",
  colors = c("#cdd5ae", "white", "#e58994"),
  hc.order = TRUE,
  hc.method = "complete",
  lab = TRUE,
  title = "",
  legend.title = "Average Correlation Coefficient",
  lab_size = 2.5
) +  
  theme(
    #legend.position = "right",         
    legend.title = element_text(size = 7),  
    axis.text.x =element_text(color="black",size=7,family="Arial" ) ,
    axis.text.y =element_text(color="black",size=7,family="Arial" ) 
  )                

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/count.cor/pa.count.cor(all.mean)","2.pptx"), width = 7, height = 7)

print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/count.cor/pa.count.cor(all.mean)","2.png"),width=7,height=7) 
  
###（3.2）Commonly identified pA sites between each pair of tools---------------------------------

res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/pair.jaccard.index(based.3seq).csv")
cor.res=read.csv("G:/scAPA/summary/temp/fig2/bech_pacount/pair.count.cor(based.3seq).csv")

index="cor"

 if(index=="jaccard") {
   data=res[res$type=="common",] %>%
     #group_by(species,protocol,method, compare.method) %>%
     group_by(protocol,method, compare.method) %>%
     summarise(
       mean_jaccard = mean(jaccard.index,na.rm = TRUE),
       num_samples = n(),
       .groups = 'drop'
     )
     
 }else if(index=="cor"){
   data=cor.res[cor.res$padj<0.05,] %>%
     #group_by(species,protocol,method, compare.method) %>%
     group_by(protocol,method, compare.method) %>%
     summarise(
       mean_cor = mean(cor,na.rm = TRUE),
       num_samples = n(),
       .groups = 'drop'
     )
 
 }
  
data$protocol<-factor(data$protocol,
                      level=c("10x" ,"CEL-seq","Drop-seq","Microwell-seq"),
                      labels=c("10X" ,"CEL-seq","Drop-seq","Microwell-seq"))



seq.col=c("10X"="#a5cc26","CEL-seq"="#A46ACD","Drop-seq"="#E78D35","Microwell-seq"="#002fa7")

data$method=factor(data$method,
                   levels = c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","MAAPER","SCAPTURE","Infernape","scUTRquant","scraps"))


##(3.2.1)line graphs--
tool.lists=unique(data$method)

plot.lists=list()
for(t in 1: length(tool.lists)){
  tool=tool.lists[t]
  
  temp.data=data[data$method==tool,]
  test = temp.data %>%
    dplyr::group_by(method,compare.method) %>%
    dplyr::summarise(
      all.mean.cor = mean(mean_cor, na.rm = TRUE),
      #sd_value = sd(value, na.rm = TRUE),
      .groups = "drop"  
    ) %>%
    arrange(desc(all.mean.cor)) 
  
  temp.data$compare.method<- factor(temp.data$compare.method,
                       levels=test$compare.method)
  
  
  p=ggplot(temp.data, aes(x = compare.method,y = mean_cor, color =protocol,fill= protocol,group = protocol)) +  
    #geom_line(linetype="dashed",size=1) +
    geom_line() +
    geom_point(size=2,shape=21) +  
    scale_color_manual(values=seq.col)+
    scale_fill_manual(values=seq.col)+
    labs(x = "", y = "correlation coefficient",
         title = tool)+
    p_theme+
    guides(fill="none",
           color="none")+
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45,hjust=1),
          # 主标题居中设置
          plot.title = element_text(
            hjust = 0.5,  
            size = 10,    
            face = "bold"  
          ),
          legend.text = element_text(size = 7))+ 
    scale_y_continuous(limits = c(0.2, 1.0),   
                       breaks = seq(0.2, 1.0, 0.2)) 
    
    
    
plot.lists[[tool]]=p
}


plot.lists <- lapply(seq_along(plot.lists), function(i) {
  plot.lists[[i]] + 
    scale_y_continuous(limits = c(0.2, 1.0)) +
    theme(
      axis.title.y = if (i %% 2 == 1) element_text() else element_blank()
    )
})



p=cowplot::plot_grid(plotlist =plot.lists,  
                   ncol = 2,
                   align = "hv",  
                   axis ="tblr" ) 



topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/count.cor/pair.count.",index,".pptx"), width = 7, height = 10)

print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig2/bech_pacount/count.cor/pair.count.",index,".png"),width=7,height=10)



  
