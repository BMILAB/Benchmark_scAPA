
library(movAPA)
library(dplyr)
library(tidyr)
library(patchwork)
library(SummarizedExperiment)
library(ggplot2)
library(eoffice)
#library(extrafont)
library(export)
library(rtracklayer)
library(GenomicFeatures)
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
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


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

sample.col=c("Mouse_Sperm"="#775360","Mouse_T_Cell"="#e8b9c5","Mouse_TIP"="#a85559","Mouse_Stem_Cell"="#976E30","Mouse_GFP1"="#ecdfa6","Mouse_GFP2"="#ecdfa6",
             "Mouse_Intestine1"="#8080A3","Mouse_Intestine2"="#8080A3","Mouse_Intestine3"="#8080A3","Mouse_Intestine4"="#8080A3","Mouse_Intestine5"="#8080A3","Mouse_Intestine6"="#8080A3",
             "Mouse_Hypothalamus1"="#bfd3e1","Mouse_Hypothalamus2"="#bfd3e1","Mouse_Bone1"="#3d85b9","Mouse_Bone2"="#3d85b9",
             "Human_Pbmc1"="#133f7f","Human_Pbmc2"="#133f7f","Human_Pbmc3"="#133f7f","Human_Pbmc4"="#133f7f","Human_Covid_Pbmc1"="#CF6862","Human_Covid_Pbmc2"="#CF6862",
             "Arabidopsis1"="#577149","Arabidopsis2"="#577149","Arabidopsis3"="#577149")

tool.col=c("scAPAtrap"="#064D4B","Sierra"="#fbd26a","SCAPE"="#ABC8E5","scAPA"="#6666CC","polyApipe"="#e8b9c5","MAAPER"="#64d9b8",
           "SCAPTURE"="#b5182b","Infernape"="#775360","scUTRquant"="#7fb80e","scraps"="#FF6B00")


heatmap.col=c("#313772","#2c4ca0","#326db6","#478ecc","#75b5dc","#fee3ce","#eabaa1","#dc917b","#d16d5b","#c44438","#b7282e")

###Create the dataset required for DTU evaluation-------------------------------------------------
##（1）Cell type annotation information----------------------------------------

anno.lists=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/anno.lists.rds")



#（2）PAC dataset----------------------------------
#(2.1)Using the conservative APA gene-----------------------------------------------
  samples=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
            "arab")
  tool.lists2=c("scAPAtrap","Sierra","SCAPE","scAPA","polyApipe","SCAPTURE","Infernape","scUTRquant","scraps")
  

  
  
  apagene=data.frame()
  apagene.lists=list()
  paclists=list()
  
  for(i in 1: 3){
    sample=samples[i]
    
    
    ##scPACdataset
    temp.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
    
    GSM=names(temp.pac)
    
    
     for(j in 1: length(tool.lists2)) {
       tool=tool.lists2[j]
       pac=temp.pac[[GSM]][[tool]]
      
        #pac=subsetPACds(pac,totPACtag = ncol(pac@counts)/20,verbose = T)
       cell.anno=anno.lists[[GSM]]
       
       ##去除barcode后面的[1-9]
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
       
       
       anno<- pac@anno %>%
         group_by(gene) %>%
         filter(n() > 1) %>%  
         ungroup()
       
       
         
         temp=data.frame(sample=GSM,
                         method=tool,
                         pac.num=nrow(pac@anno),
                         APA.gene.num=length(unique(anno$gene)),
                         APA.gene=unique(anno$gene))
        
      
     
       apagene=rbind(apagene,temp)
       apagene.lists[[GSM]]=apagene
     
    cat(paste0(GSM,": ",j,"/",length(tool.lists2)))
     }
  }
  
  #write.csv(sample.data,file = paste0("G:/scAPA/summary/temp/fig2/","three.scapa.number.csv"),row.names = FALSE)
  saveRDS(apagene.lists,file = paste0("G:/scAPA/summary/temp/fig3/","apagene.lists.rds"))
  saveRDS(paclists,file = paste0("G:/scAPA/summary/temp/fig3/","sc.paclists.lists.rds"))
  
  ##Filtering conserved APA genes------------------------- 

  com.gene=list()
  
  for(i in 1: length(apagene.lists)){
    sample=names(apagene.lists)[i]
    
    bulk.pac=bulk.paclists[[i]]
    bulk.anno<- bulk.pac@anno %>%
      group_by(gene) %>%
      filter(n() > 1) %>%  
      ungroup()
    
    info=apagene.lists[[sample]] %>%
      dplyr::count(`APA.gene`,sort=TRUE)
    ##keeping APA genes identified by at least six different tools.
    info=subset(info,info$n >=6)
    
    com=intersect(info$APA.gene,bulk.anno$gene)
    
    com.gene[[sample]]=com
    
  }
  
  saveRDS(com.gene,file ="G:/scAPA/summary/temp/fig3/com.gene.rds")  
  
  
  
##2.Filter pA based on expression levels, retain APA genes, and then perform DTU analysis.---------------------
  
  sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
  #GSMlists=unique(sample.data$Sample_ID)
  anno.lists=readRDS("G:/scAPA/summary/temp/fig2/bech_pacount/anno.lists.rds")
  
  
  GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP",
             "3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490691")
 
  
  
  anno.lists=subset(anno.lists, names(anno.lists) %in% GSMlists)
  
  res.lists=list()
  data.res=data.frame()
  
  for(i in c(7: 14)){
    GSM <-  GSMlists[i]
    
    species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
    protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
    sample<- unique(sample.data[sample.data$Sample_ID==GSM,]$Sample)
    
    if (i == 1) {
      #celltypes=c("SC","RS","ES")
      groups=c("ES","SC")
      
    } else if (i == 2) {
      #celltypes=c("naive_Tcell","activated_Tcell")
      groups=c("activated_Tcell","naive_Tcell")
      
    } else if (i == 3) {
      #celltypes=c("Sham_Fibroblasts","Sham_leukocytes","Sham_EC","MI_leukocytes")
      groups=c("Sham_Fibroblasts","Sham_EC")
      
    } else if (i == 4) {
     
      groups=c("ESC","MEF")
      
    } else if (i == 5 |i == 6 ) {
     
      groups=c("F-Act","F-SH/F-SL")
      
    } else if (i %in% c(7,8,9,10,11) ) {
      #CD14+_monocyte +CD16+_monocyte =Monocyte
      groups=c("Monocyte","B_cells")
    } else if (i == 12) {
      #CD4+ T  + CD8+ T= T_cells
      groups=c("Monocyte","T_cells")
      
    } else if (i %in% c(13,14)) {
     
      groups=c("Atrichoblast","Trichoblast")
      
    } else {
      stop("Error：差异分析的样本数为14，i 在14以内")
      
    }
    
    
    ##scPACdataset
    temp.pac=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
    
    #GSM=names(temp.pac)
    
    #Cell type annotation information
   
      
      for (i in c(5,6)) {
        anno.lists[[i]]$CellType <- ifelse(
          anno.lists[[i]]$CellType %in% c("F-SH", "F-SL"), 
          "F-SH/F-SL",                                   
          anno.lists[[i]]$CellType                       
        )
      }
      
   
      
      for(i in c(7,8,9,10)){
        anno.lists[[i]]$CellType <- gsub("CD14\\+_monocyte", "Monocyte", anno.lists[[i]]$CellType)
        anno.lists[[i]]$CellType <- gsub("CD16\\+_monocyte", "Monocyte", anno.lists[[i]]$CellType)
      }  
      
   
     
      for(i in c(11,12)){
        anno.lists[[i]]$CellType <- gsub("B cell", "B_cells", anno.lists[[i]]$CellType)
        anno.lists[[i]]$CellType <- gsub("CD14Mono", "Monocyte", anno.lists[[i]]$CellType)
        anno.lists[[i]]$CellType <- gsub("CD16Mono", "Monocyte", anno.lists[[i]]$CellType)
        anno.lists[[i]]$CellType <- gsub("CD4\\+ T", "T_cells", anno.lists[[i]]$CellType)
        anno.lists[[i]]$CellType <- gsub("CD8\\+ T", "T_cells", anno.lists[[i]]$CellType)
        
      } 
    
      
    tool.lists=names(temp.pac[[GSM]])
    tool.lists=setdiff(tool.lists, "MAAPER")
    
    for(j in 1: length(tool.lists)) {
      tool=tool.lists[j]
      pac=temp.pac[[GSM]][[tool]]
      
      #pac=subsetPACds(pac,totPACtag = ncol(pac@counts)/20,verbose = T)
      cell.anno=anno.lists[[GSM]]
      ##(1)Add cell type annotation information---------------------------------
      
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
      pac@colData$celltype=cell.anno$CellType[loc]
      
      if(anyNA(loc)==TRUE){
        stop(paste0(GSM," (",tool,")",": 细胞类型注释信息不匹配"))
      }
      
      
      ##(2)Filter out pAs with total expression< 10 in thesample group and retain APA genes---------------------
      sub.pac=subsetPACds(pac,group="celltype",conds = groups, totPACtag=10,verbose = T)
      
     
      #paclists[[GSM]][[tool]]=sub.pac
      
      
     # anno<- sub.pac@anno %>%
     #   group_by(gene) %>%
     #   filter(n() > 1) %>%  # 保留转录本数 >1 的基因
     #   ungroup()
      
      
      #temp=data.frame(sample=GSM,
      #                method=tool,
     #                 pac.num=nrow( sub.pac@anno),
     #                 APA.gene.num=length(unique(anno$gene)),
     #                 APA.gene=unique(anno$gene))
      
      
      
     # apagene=rbind(apagene,temp)
     # apagene.lists[[GSM]]=apagene
      
      ##（3）DE APA Analysis (Based on constructed pseudobulk subgroups)---------------------------
      ##DEXseq requires that gene names or feature names do not contain special characters such as “:” or spaces
      row.names(sub.pac@anno)=paste0(sub.pac@anno$chr,";",sub.pac@anno$strand,";",sub.pac@anno$coord)
      row.names(sub.pac@counts)=row.names(sub.pac@anno)
      
      
      counts=sub.pac@counts
      
      tx2gene=data.frame(tx_id=row.names(sub.pac@anno),
                         gene_id=sub.pac@anno$gene
      )
      
      samples=data.frame(sample_id=sub.pac@colData$barcode,
                         condition=sub.pac@colData$celltype)
      
      cts <- tibble::rownames_to_column(as.data.frame(counts), var = "tx_id")
      cts <- merge(cts, tx2gene, by = "tx_id") %>%
        dplyr::select(gene_id, tx_id, everything())
      
      #过滤掉只对应唯一feature_id的gene_id
      cts <- cts %>%
        group_by(gene_id) %>%
        filter(n_distinct(tx_id) > 1) %>% 
        ungroup() %>%
        distinct()
      
      cts <- as.data.frame(cts)
      tx2gene <- cts[, c(2, 1)]
      
      #constructed pseudobulk subgroups---------------------
      meta.data=data.frame(row.names = samples$sample_id,
                           condition=samples$condition)
      row.names(cts)=cts$tx_id
      cts=cts[,-c(1,2)]
      
      pseudo.bulk.cts=create_pseudo_bulk(counts=cts, samples=meta.data,
                                         population.1 = groups[1], 
                                         population.2 = groups[2], 
                                         num.splits = 6, 
                                         seed.use = 1,
                                         replicates.1 = NULL,
                                         replicates.2 = NULL)
      
      
      pseudo.bulk.cts <- tibble::rownames_to_column(as.data.frame(pseudo.bulk.cts), var = "tx_id")
      pseudo.bulk.cts <- merge(pseudo.bulk.cts, tx2gene, by = "tx_id") %>%
        dplyr::select(gene_id, tx_id, everything())
      
      samples=data.frame(sample_id=colnames(pseudo.bulk.cts[,-c(1,2)]),
                         condition=colnames(pseudo.bulk.cts[,-c(1,2)]))
      
      samples$condition <- ifelse(grepl("^Population1_", samples$condition), groups[1],
                                  ifelse(grepl("^Population2_",samples$condition), groups[2],
                                         samples$condition))
      
      tx2gene <- pseudo.bulk.cts[, c(2, 1)]
      
      
      input_data <- PreInput(counts=pseudo.bulk.cts, samples, tx2gene)
      
      
      ###DTU Analysis------------------------------------------
      
      temp.res <- tryCatch({
        run_DEXSeq(input_data,conditions=groups)
      }, 
      error = function(e) {
        message(paste0("Error in ", GSM, "(", tool, "):  ", e$message))
        return(NULL)  
      })
      
     
      if (is.null(temp.res)) {
        next
      }else{
        
        temp.res$log2FC=temp.res[,10]
        
        temp.res$contr=paste0(groups[1],"-",groups[2])
        
        temp.res <- temp.res %>%
          mutate(change = case_when(
            is.na(log2FC) | is.na(padj) ~ "not",
         
            (abs(log2FC) > 0.585 & padj < 0.05) ~ "sig",
            
            TRUE ~ "not"
          )) 
      
      
      
      filter.res=subset(temp.res,temp.res$change=="sig")
      
      res.lists[["all"]][[GSM]][[tool]]=temp.res
      res.lists[["sig"]][[GSM]][[tool]]=filter.res
      
      
      
      temp=data.frame(species=species,
                      protocol=protocol,
                      sample=sample,
                      sample_id=GSM,
                      tool=tool,
                      gene.num=length(unique(temp.res$gene_id)),
                      apa.num=nrow(temp.res),
                      sig.gene.num=length(unique(filter.res$gene_id)),
                      sig.apa.num=nrow(filter.res),
                      
                      contr=paste0( groups[1], "-", groups[2], sep = ""))
      
      
      data.res=rbind(data.res,temp)
      }
     
      message(paste0(GSM,"(",tool,")",": ",j,"/",length(tool.lists)))
    }
  }
  
  #write.csv(sample.data,file = paste0("G:/scAPA/summary/temp/fig2/","three.scapa.number.csv"),row.names = FALSE)
  write.csv(data.res,file="G:/scAPA/summary/temp/fig3/deXseq/dtu.sta_pseudo.bulk.csv",row.names = FALSE)
  saveRDS(res.lists,file="G:/scAPA/summary/temp/fig3/dexseq/dtu.res_pseudo.bulk.rds")  
  
  #res.lists=lapply(res.lists[c("all", "sig")], function(x) x[1:4])

 
  
  ##3.Integrate DEAPA genes from different tools to construct a core consensus gene set-------------
  
 ##(3.1)Compare genes with differential APA usage--------------------------------------------
  
   anno.deapa=readRDS("G:/scAPA/summary/temp/fig3/dexseq/dtu.res_pseudo.bulk.rds")
  
  #Retain samples containing a higher number of DEAPA transcripts
  #GSMlists=names(anno.deapa[["sig"]])
 GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
                       "4k","5k","GSM4712885","GSM4712907")
 
 ##Criteria for filtering consensus genes
 
 min.wf=1.6
 min.tool=4
 
  
  com.de.genelists=list()
  deapa.gene.lists=list()
  for( i in 1 :8){
    GSM=GSMlists[[i]]
    #sample=samples[i]
    
    
    tool.lists=names(anno.deapa[["sig"]][[GSM]])
    
   
    for(j in 1: length(tool.lists)){
      
      tool=tool.lists[j]
      deapa=anno.deapa[["sig"]][[GSM]][[j]]
      
      ##（1）For each gene, retain the padj corresponding to the pA with padj < 0.05 and the largest |log2FC|
      deapa <- deapa %>%
        group_by(gene_id) %>%               
        arrange(desc(abs(log2FC)), padj,  .by_group = TRUE) %>%  
        slice_head(n = 1) %>%                                        
        ungroup()                           
      
      ##（2）Calculate significance rankings
      deapa <- deapa %>%
        arrange(padj, desc(abs(log2FC))) %>%   
        mutate(padj.rank = row_number())   
      deapa$PR=(deapa$padj.rank/nrow(deapa))* 100 
      
      deapa=deapa[,c("gene_id","log2FC","padj","padj.rank","PR")]
      
      
      deapa.gene.lists[[GSM]][[tool]]=deapa
      
    }
    
    ##（3)Weighted Frequency, WF
    # WF_gene = Σ (1 - PR) × I(DEG)
    wf.res=calculate_wf(deapa.gene.lists[[GSM]])
    
    
    ##(4)Filter consensus genes based on WF and supporting_tools numbers
    
    temp=as.data.frame(table(unlist(lapply(deapa.gene.lists[[GSM]], function(df) df[, 1]))))
    colnames(temp)=c("gene","tool.num")
    
   
    res=left_join(wf.res,temp, by = "gene")
    
    #查看不同分位数所对应的wf值
    #quantile(res$wf, probs = c(0.6, 0.8)) 
    
    sub.res=subset(res,res$wf>=min.wf & tool.num>=min.tool)
    
    
    com.de.genelists[[GSM]]=sub.res
    
    message(paste0(GSM,": ",i,"/",length(GSMlists)))
    
  }
  
 
  saveRDS(deapa.gene.lists,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/old/deapa.gene.lists.rds")  
  saveRDS(com.de.genelists,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/old/com.de.genelists.rds")  
  
######################################################################################
 
  
  ##（3.2）Compare genes with differential APA usage and, 
  ##if possible, ensure that the direction of regulation is consistent.-------------------------------- 
  
  dtu=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/dtu.res_pseudo.bulk.rds")
  
  ##Retain samples containing a higher number of DEAPA transcripts
  #GSMlists=names(anno.deapa[["sig"]])
  GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
             "4k","5k","GSM4712885","GSM4712907")
  
  ##Criteria for filtering consensus genes
  
  min.wf=1.6
  min.tool=4
  
  com.de.genelists=list()
  deapa.gene.lists=list()
  for( i in 1 :8){
    GSM=GSMlists[[i]]
    
    
    tool.lists=names(dtu[["all"]][[GSM]])
    
    for(j in 1: length(tool.lists)){
      tool=tool.lists[j]
      
      deapa=dtu[["all"]][[GSM]][[tool]]
      ###(1)To determine the direction of regulation for pA sites with the largest |log2FC| value 
      ##and a p-adjusted (padj) value less than 0.05, follow these steps:
      
      
      deapa= deapa[!is.na(deapa$log2FC) | !is.na(deapa$padj), ]
      
      
      deapa <- deapa %>%
        separate(tx_id, into = c("chr", "strand", "coord"), sep = ";", remove = FALSE)%>%
        mutate(
          coord = as.numeric(coord))
      
      deapa <- deapa %>%
        group_by(gene_id) %>%
        mutate(
          DEAPA.gene = ifelse(any(change == "sig"), TRUE, FALSE),
          rank = ifelse(
            strand == "+",  
            rank(coord, ties.method = "first"),  
            rank(-coord, ties.method = "first")  
          ),
          label = paste0(gene_id, "_", rank)  
          
        ) %>%
        ungroup() %>%
        arrange(gene_id, rank)
      
    
      deapa=subset(deapa,deapa$DEAPA.gene==TRUE)
      
      if(nrow(deapa)==0){
        
        col <- c("gene_id","gene","log2FC","padj","padj.rank","PR")
        
        df <- data.frame(matrix(ncol = length(col), nrow = 0))
        colnames(df) <- col
        
        deapa.gene.lists[[GSM]][[tool]]=df
      }else{
        
        deapa <- deapa %>%
          
          filter(change == "sig")%>% 
          
          group_by(gene_id) %>%
          mutate(
            max_abs_estimate_row = which.max(abs(log2FC)),
           
            top_rank = rank[max_abs_estimate_row],
            
            top_estimate_sign = sign(log2FC[max_abs_estimate_row]),
            
            DEAPA.gene.type = case_when(
              top_rank == 1 & top_estimate_sign > 0 ~ "shorted",
              top_rank == 1 & top_estimate_sign < 0 ~ "lengthened",
              top_rank != 1 & top_estimate_sign > 0 ~ "lengthened",
              top_rank != 1 & top_estimate_sign < 0 ~ "shorted",
              TRUE ~ NA  
            )
          ) %>%
          dplyr::select(-max_abs_estimate_row, -top_rank, -top_estimate_sign) %>%  # 移除临时列
          ungroup()
        
        deapa$res=paste0(deapa$gene_id,"_",deapa$DEAPA.gene.type)
        
        ####（2）For each gene, retain the padj corresponding to the pA with padj < 0.05 and the largest |log2FC|
        deapa <- deapa %>%
          group_by(gene_id) %>%               
          arrange(desc(abs(log2FC)), padj,  .by_group = TRUE) %>%  
          slice_head(n = 1) %>%                                         
          ungroup()                          
        
        ##（3）Calculate Significance Rankings
        deapa <- deapa %>%
          arrange(padj, desc(abs(log2FC))) %>%  
          mutate(padj.rank = row_number())   
        deapa$PR=(deapa$padj.rank/nrow(deapa))* 100 
        
        deapa=deapa[,c("res","gene_id","log2FC","padj","padj.rank","PR")]
        colnames(deapa)[1:2]=c("gene_id","gene")
        
        deapa.gene.lists[[GSM]][[tool]]=deapa
        
      }
    }
    
    ##（3)Weighted Frequency, WF
    # WF_gene = Σ (1 - PR) × I(DEG)
    wf.res=calculate_wf(deapa.gene.lists[[GSM]])
    
    
    ##(4)Screen consensus genes based on wf (weighted frequency) and the number of tools supporting them.
    
    
    temp=as.data.frame(table(unlist(lapply(deapa.gene.lists[[GSM]], function(df) df[, 1]))))
    colnames(temp)=c("gene","tool.num")
    
    
    res=left_join(wf.res,temp, by = "gene")
    
    sub.res=subset(res,res$wf>=min.wf & tool.num>=min.tool)
    
    
    com.de.genelists[[GSM]]=sub.res
    
    message(paste0(GSM,": ",i,"/",length(GSMlists)))
    
  }
  
  saveRDS(deapa.gene.lists,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.gene.lists.rds")  
  saveRDS(com.de.genelists,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/com.de.genelists.rds")   
  
  
  
  
  
 
  
###4.Count the number of overlapping DEAPA genes--------------------------------

  deapa.gene.lists=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.gene.lists.rds")
  com.de.genelists=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/com.de.genelists.rds")
  
 
  
  
  GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
             "4k","5k","GSM4712885","GSM4712907")
  
  
  res=data.frame()
  for( i in 1 :8){
    GSM=GSMlists[[i]]
   
    
    tool.lists=names(deapa.gene.lists[[GSM]])
    
    for(j in 1: length(tool.lists)){
      tool=tool.lists[j]
      
      deapa=deapa.gene.lists[[GSM]][[tool]]
      
     com.gene=com.de.genelists[[GSM]]
       
        temp.res=data.frame(sample=GSM,
                            tool=tool,
                            com.gene.num=nrow(com.gene),
                            de.gene.num=nrow(deapa),
                            overlap.gene.num=length(intersect(deapa$gene_id,com.gene$gene))
                      
                            
        )
        
        ##Calculate precision and recall
        temp.res <-temp.res %>%
          group_by(sample, tool) %>%
          mutate(
            precision = overlap.gene.num / de.gene.num,
          
            recall = overlap.gene.num / com.gene.num
          ) %>%
          #fill(precision, recall, .direction = "downup") %>%
          ungroup()
        
        
        res=rbind(res,temp.res)
        
        message(paste0(GSM,"(",tool,")",": ",j,"/",length(tool.lists)))
        
      }
      
    }
  

 # write.csv(res,file="G:/scAPA/summary/temp/fig3/dexseq/fil.by.exp/deapa.res.csv",row.names = FALSE)  
 # write.csv(res,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/old/deapa.res.csv",row.names = FALSE)
  write.csv(res,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.res.csv",row.names = FALSE)
 
  
###4.plot-----------------------------------------------------------------------
###（1）Number of validated DEAPA genes--------------------------------------------------
  res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.res.csv")
  
  head(res)
 
  sub.res =res[,c(1,2,4,5)] %>% 
    mutate(uni.de.gene.num=de.gene.num-overlap.gene.num)
    
  input.data<- sub.res%>%
    pivot_longer(
      cols = c(overlap.gene.num, uni.de.gene.num),  
      names_to = "type",               
      values_to = "value"              
    ) %>%
    mutate(
      type = case_when(
        type == "overlap.gene.num" ~ "overlap",     
        type == "uni.de.gene.num" ~ "unique" 
      )
    )
  
  
  input.data$sample <- factor(input.data$sample,
                                 levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
                                          "4k","5k","GSM4712885","GSM4712907"),
                                 labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell",
                                            "Human_Pbmc2","Human_Pbmc3","Human_Covid_Pbmc1","Human_Covid_Pbmc2"))
                              
  
  ##tool order--------------
  
  
  #order=res %>% 
  #  group_by(tool)%>%
  #  summarise(mean.precision=mean(precision, na.rm = TRUE),
  #         mean.recall=mean(recall, na.rm = TRUE))%>%
  #  arrange(desc(mean.recall))
  
  #order=res %>% 
  #  group_by(tool)%>%
  #  summarise(mean.de.gene.num=mean(de.gene.num, na.rm = TRUE))%>%
  #  arrange(desc(mean.de.gene.num))
  
  order <- res %>%
    group_by(sample) %>%
    mutate(rank = rank(-de.gene.num, ties.method = "min")) %>%
    ungroup() %>%
    group_by(tool) %>%
    summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
    arrange(mean_rank)
  

  
  
  
  input.data$tool=factor(input.data$tool,
                             levels = rev(order$tool))
  
  input.data$type <-factor(input.data$type,
                               levels =rev(c("overlap","unique")),
                               labels =rev( c("Overlapping","Unique"))
  )
  
 
  
  cols=c("Unique"="#b1cde5","Overlapping"="#004f9b") 
  
  cols=c("Unique"="#bedeeb","Overlapping"="#437c93") 
 
  ##柱状图种文本标签所显示的结果
  text.data=subset(input.data,input.data$type=="Overlapping")
  text.data$y.label=text.data$value *0.5
  
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
  
  
  p=ggplot(data=input.data,
           aes(x=tool,  y=value, fill=type)) +
    geom_bar(stat="identity",
             width = 0.7,                 
             position = "stack" 
    )+
    
    geom_text(data=text.data,aes(y=y.label, label=paste0(value)), 
              #position = position_stack(vjust = 0.5),  
              #vjust = -0.2,
              #hjust = -0.1,
              color="#fed6a2", #"#fed6a2"
              size=2)+
    p_theme+
    labs(x=NULL,y="No. of DE APA Gene",title = "")+
    guides(fill=guide_legend(title=NULL))+
    scale_fill_manual(values=cols)+
    
    theme(#legend.position = "top",
      #axis.text.x = element_text(angle = 45,hjust=1),
      #strip =FALSE, 
      strip.text.x = element_text(size=8),
      element_rect(fill = "grey"),
      panel.spacing.x = unit(0.5, "cm"), 
      )+
    coord_flip()+
    facet_wrap(~ sample,scales="free_x",nrow=2)
    #facet_wrap(~ sample,nrow=2)
  
  
  topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/deapa.num2.pptx"), width = 8, height = 6)
 
  #保存为png
  print(p)
  graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/deapa.num2.png") ,width=8,height=6)
  
  
  
 ##(2)Precision and Recall-------------------------------------------------------
  index="Precision"
  
  
  res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/deapa.res.csv")
  res$precision[is.na(res$precision)] <- 0
  
  input.data =res %>% 
    group_by(tool)%>%
    mutate(mean.precision=mean(precision, na.rm = TRUE),
           mean.recall=mean(recall, na.rm = TRUE))%>%
    arrange(desc(mean.precision))
  
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
  
  
  
  input.data$tool=factor(input.data$tool,
                         levels = rev(order$tool))
  
  input.data$sample <- factor(input.data$sample,
                              levels=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
                                       "4k","5k","GSM4712885","GSM4712907"),
                              labels = c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell",
                                         "Human_Pbmc2","Human_Pbmc3","Human_Covid_Pbmc1","Human_Covid_Pbmc2"))
  
  
###Lollipop plot----------------------------------------------------------------  
  p=ggplot(input.data, aes(x =tool, y = recall,fill=tool,group=tool)) +
    geom_segment(
      aes(x = tool, xend = tool, y = 0, yend = recall),
      color = "grey50", 
      linewidth = 1      
    ) +
    # 绘制圆点
    geom_point(
      aes(size = overlap.gene.num),  
      #size=5,
      color = "grey50",      
      #fill = "white",        
      shape = 21,             
      stroke = 1.5,            
      alpha=0.7
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
  p_theme+
    labs(x=NULL,y=index,title = "")+
    guides(size=guide_legend(title="Overlapping DE APA Gene"))+ 
    guides(fill='none')+
    scale_fill_manual(values=tool.col)+
    
    theme(legend.position = "top",
      strip.text.x = element_text(size=8),
      element_rect(fill = "grey"),
      panel.spacing.x = unit(0.5, "cm"),  
     
    )+
    coord_flip()+
    facet_wrap(~ sample,nrow=2)
  
  topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/",index,".pptx"), width = 8, height = 6)
  
  #保存为png
  print(p)
  graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/",index,".png") ,width=8,height=6)
  
###violin plot----------------------------------------------------------- 
  input.data$tool=factor(input.data$tool,
                         levels = order$tool)
  
 p= ggplot(input.data, aes(x =tool, y = precision,group=tool)) +
    geom_violin( aes(color = tool),  
                 fill = NA,         
                 trim = FALSE,
                 linewidth = 1   ) + 
    
    geom_point(
      #aes(size = overlap.gene.num),  
      size=1.5,
      color = "grey50",      
      #fill = "grey50",         
      shape = 21,            
      stroke = 1.2,            
      alpha=0.7
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
    p_theme+
    labs(x=NULL,y=index,title = "")+
    guides(size=guide_legend(title="Overlapping DE APA Gene"))+ 
    guides(fill='none',color='none')+
    scale_color_manual(values=tool.col)+
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 45,hjust=1),
          strip.text.x = element_text(size=7),
          element_rect(fill = "grey"),
          panel.spacing.x = unit(0.5, "cm"),  
    )
   
  
  topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/",index,"(violin.plot).pptx"), width = 5, height = 3)
  
  print(p)
  graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/",index,"(violin.plot).png") ,width=5,height=3) 
  
 
  
## (3)The consistency of DEAPA gene identification among different methods------------------------------------ 
  anno.deapa=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.gene.lists.rds")
  
  compare.method="padj" #padj/log2FC
  
  #GSMlists=names(anno.deapa)
  GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
             "4k","5k","GSM4712885","GSM4712907")
  
  res=data.frame()
  res.jaccard=data.frame()
  for( i in 1 :8){
    GSM=GSMlists[[i]]
    #sample=samples[i]
    
    tool.lists=names(anno.deapa[[GSM]])
    #tool.lists=names(anno.deapa[["sig"]][[GSM]])
    
    for(j in 1: length(tool.lists)){
      
      tool=tool.lists[j]
      
      deapa=anno.deapa[[GSM]][[j]]
      
  
    for(k in tool.lists){
      temp.deapa=anno.deapa[[GSM]][[k]]
      
     
      #Keep the commonly identified DEAPA genes
      info=intersect(deapa$gene_id, temp.deapa $gene_id)
      
      sub.deapa=as.data.frame(deapa)
      row.names(sub.deapa)=sub.deapa$gene_id
      
      sub.temp.deapa=as.data.frame(temp.deapa)
      row.names(sub.temp.deapa)=sub.temp.deapa$gene_id
      
      sub.deapa=sub.deapa[info,]
      sub.temp.deapa=sub.temp.deapa[info,]
      
      ##（3.1）Jaccard index--------------------
      all.gene=union(deapa$gene_id, temp.deapa $gene_id)
      
      temp.jaccard=data.frame(sample=GSM,
                              method=tool,
                              compare.method=k,
                              overlap.gene.num=length(info),
                              all.gene.num=length(all.gene),
                             jaccard.index= length(info)/ length(all.gene)
                              
      )
      
      res.jaccard=rbind( res.jaccard,temp.jaccard)
      
      ###（3.2）Comparison of log2FC/padj correlation------------------------------
      #If the number of commonly identified DE APA genes is less than 10, then the correlation will not be calculated
      if(length(info)<10){
        next
      }else{
        
        ###log2FC corr---------------------------   
        
        if(compare.method=="log2FC"){
          
          temp=data.frame(row.names = row.names(sub.deapa),
                          log2fc=abs(sub.deapa$log2FC),
                          temp_log2fc=abs(sub.temp.deapa$log2FC)
                          
          )
          
          
          temp<- temp %>%
            mutate(diff = log2fc-temp_log2fc) %>%
            arrange(abs(diff))  %>% 
            #slice_head(n = 600)
            slice_head(prop= 1)
          
          
          
          cor=psych::corr.test(temp$log2fc,  temp$temp_log2fc,
                               use = "pairwise", 
                               method="pearson", 
                               adjust = "fdr"   
          )
          
          
          temp.res=data.frame(sample=GSM,
                              method=tool,
                              compare.method=k,
                              cor=cor$r,
                              p.adj=cor$p.adj,
                              gene.num=cor$n
                              
          )
        
          
          #padj corr-------------------------------------------------
          
        } else if(compare.method=="padj"){
          
          sub.temp.deapa=sub.temp.deapa[row.names(sub.deapa),]
          
          
          temp=data.frame(row.names = row.names(sub.deapa),
                          padj.rank=sub.deapa$padj.rank,
                          temp_log2fc=sub.temp.deapa$padj.rank
          )
          
          cor=psych::corr.test(temp$padj.rank,  temp$temp_log2fc,
                               use = "pairwise", 
                               method="spearman", 
                               adjust = "fdr"   
          )
          
          temp.res=data.frame(sample=GSM,
                              method=tool,
                              compare.method=k,
                              cor=cor$r,
                              p.adj=cor$p.adj,
                              gene.num=cor$n
                              
          )
          
        }
        
        res=rbind(res,temp.res)
     
      }
    }
    
    message(paste0(GSM,": ",j,"/",length(tool.lists)))
  }
 }
  
  write.csv(res.jaccard,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/res.jaccard(precies.gene).csv",row.names = FALSE)
  #write.csv(res,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/tool_cor_log2fc(precies.gene).csv",row.names = FALSE)
  write.csv(res,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/tool_cor_padj(precies.gene).csv",row.names = FALSE)
  
####（4）plot---------------------------------------------------------------------------- 
##（4.1）Heatmap of the jacard index----------------------------------------------------------  
  
  res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/res.jaccard(precies.gene).csv") 
  
  
###Calculate the average
  data= res %>%
    group_by(method, compare.method) %>%
    summarise(
      mean_jaccard = mean(jaccard.index,na.rm = TRUE),
      num_samples = n(),
      .groups = 'drop'
    )
  
  
 
  cor=data[,c(1,2,3)]
  cor=tidyr::pivot_wider(cor, 
                         names_from = compare.method, 
                         values_from = mean_jaccard)

  #cor[is.na(cor)]=0
 
  
  cor=as.data.frame(cor)
  row.names(cor)=cor$method
  cor=cor[,-1]
  
  hc <- hclust(as.dist(1 - cor) )  
  order <- hc$labels[hc$order]
  
  
  data$method=factor(data$method,
                     levels = order)
  data$compare.method=factor(data$compare.method,
                             levels = order)
  ##p1:Create a Clustering Dendrogram------------------------------------------------
  
  library(dendextend)
  library(ggdendro)
  dend <- as.dendrogram(hc) %>%
    set("branches_k_color", k = 3, 
        value = c("#E6CDC6", "#a0afb6", "#B3B69C")) %>%  
    set("labels_cex", 0.8)             
  
  plot(dend)
  #cluster_groups <- cutree(hc, k = 3)  
 
  dend_data <- ggdendro::dendro_data(dend)       
  
  
  
  dend_plot <- ggplot(segment(dend_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    #geom_segment(aes(x = x, y = y, xend = xend, yend = yend,color = factor(cluster_groups))  
    #) +
    #scale_color_manual(values = c("#E6CDC6", "#a0afb6", "#B3B69C")) + 
    theme_void() +
    scale_y_reverse(expand = c(0, 0))  
  
 
  
  ###p2:Heatmap plot--------------------------------------------------------------
  
  col2=brewer.pal(9,"Blues")
  
  
  heatmap_plot=ggplot(data, aes(x = method, y =compare.method, fill =mean_jaccard)) + #fill =-1*log10(p.adj)
    geom_tile(color="white",
              lwd = 0.8, 
              linetype = 1) + 
    scale_fill_gradientn(
      colors =col2)+  
    geom_text(aes(label = round(mean_jaccard, 2)), 
              size = 3, 
              color = "black") + 
    theme_minimal()+
    theme_bw()+
    labs(x=element_blank(),y=element_blank(),
         title = "",
         fill="Average Jaccard Index",
         size="Sample Number") + 
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
  
  
  p<- heatmap_plot /  dend_plot + 
    plot_layout(heights = c(4, 1))  
  
  
  
    topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/tool.jaccard(precies.gene).pptx"), width = 6, height = 6)
    
    #保存为png
    print(p)
    graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/tool.jaccard(precies.gene).png") ,width=6,height=6)  
    

  
### (4.2)Correlation Heatmap for padj/log2fc-------------------------------------------------------------------------
  index="padj"  #log2fc
  
  #res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/tool_cor_log2fc(precies.gene).csv")
  res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/tool_cor_padj(precies.gene).csv")
 
 
  res$group=paste0(res$method,"_",res$compare.method)
  table(res$group)
  
  data<- res %>% 
    filter(gene.num >= 10, p.adj < 0.05)
  
  
  data= data %>%
    group_by(method, compare.method) %>%
    summarise(
      mean_cor = mean(cor),
      num_samples = n(),
      .groups = 'drop'
    )
  
  
  

  cor=data[,c(1,2,3)]
  cor=tidyr::pivot_wider(cor, 
              names_from = compare.method, 
              values_from = mean_cor)
 
  cor[is.na(cor)]=0
  
  
  
  cor=as.data.frame(cor)
  row.names(cor)=cor$method
  cor=cor[,-1]
  
  hc <- hclust(as.dist(1 - cor) )  
  order <- hc$labels[hc$order]
  
 # plot(hc,hang = -1, main = "Cluster Dendrogram")  # hang=-1让所有标签对齐
  
  
  data$method=factor(data$method,
                     levels = order)
  data$compare.method=factor(data$compare.method,
                     levels = order)
  
 
  col <- colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
  col2<-colorRampPalette(c("#663d74","white","#fac03d"), alpha = TRUE)
  
  ##p1:Create a clustering dendrogram----------------------------------------------------------
  
  library(dendextend)
  library(ggdendro)
  dend <- as.dendrogram(hc) %>%
    set("branches_k_color", k = 3, 
        value = c("#E6CDC6", "#a0afb6", "#B3B69C")) %>%  
    set("labels_cex", 0.8)              
  
  
  

  dend_data <- ggdendro::dendro_data(dend)      
  

  
  dend_plot <- ggplot(segment(dend_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    theme_void() +
    scale_y_reverse(expand = c(0, 0)) 
  
  
  
  ####p2:heatmap plot---------------------------------------------------------------
  
  input.data=cor
  input.data$method=row.names(input.data)
  input.data=input.data%>%
    pivot_longer(
      cols = -method,  
      names_to = "compare.method",  
      values_to = "mean_cor"  
    ) 
  input.data$mean_cor[input.data$mean_cor == 0] <- NA
  
  
  
  
  input.data$group=paste0(input.data$method,"_",input.data$compare.method)
  data$group=paste0(data$method,"_",data$compare.method)
  
  loc=match(input.data$group,data$group)
  input.data$num_samples=data$num_samples[loc]
  
  #指定工具间的顺序
  input.data$method=factor(input.data$method,
                     levels = order)
  input.data$compare.method=factor(input.data$compare.method,
                             levels = order)
  input.data$num_samples[is.na(input.data$num_samples)] <- 1
  
  
  input.data <- input.data %>%
    mutate(
      is_significant = ifelse(is.na(mean_cor), "Not significant", "Significant") 
    )
  
  
  heatmap_plot=ggplot(input.data, aes(x = method, y =compare.method, fill =mean_cor)) + #fill =-1*log10(p.adj)
    geom_tile(aes(fill = rectheat),color="black",fill="white") + 
    geom_point(
      aes(size = num_samples,        
          fill =mean_cor,          
          shape = is_significant ), 
      color = "#909fac",     
      stroke = 0.8         
    ) +
   
    scale_size_continuous(
      range = c(5, 10) ,   
      guide = guide_legend(         
        override.aes = list(shape = 22))
      #,guide = "none"         
    ) + 
    
    scale_shape_manual(
      values = c("Not significant" = 8, "Significant" = 22),
      guide = guide_legend(title = "")
    ) +
    scale_fill_gradientn(
      colors =col(8))+  #col2(10)
    geom_text(aes(label =ifelse(is_significant=="Not significant", "*", round(mean_cor, 2))), 
             size = 3, 
              color = "black") + 
    theme_minimal()+
    theme_bw()+
    labs(x=element_blank(),y=element_blank(),
         title = "",
         fill="Average Correlation Coefficient",
         size="Sample Number") + 
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
  
  
  p<- heatmap_plot /  dend_plot + 
    plot_layout(heights = c(4, 1))  
  
  
  topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/tool.",index,".cor(precies.gene)2.pptx"), width = 7, height = 6)
  
  #保存为png
  print(p)
  graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/tool.",index,".cor(precies.gene)2.png") ,width=7,height=6)  
  
 
###(5)GO Analysis------------------------------------------------------- 
  
  
  ##构建gene_id与gene_ensembl的对应关系文件
  ref.genelists=list()
  reference.file <- "G:/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf"
  ref.gene<-parseGenomeAnnotation(reference.file)
  ref.genelists[["Mouse"]]=ref.gene[["anno.frame"]]
  
  reference.file <- "G:/ref_file/GRch38fa+gtf/v34/gencode.v34.annotation.gtf"
  ref.gene<-parseGenomeAnnotation(reference.file)
  ref.genelists[["Human"]]=ref.gene[["anno.frame"]]
  

  deapa.gene.lists=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.gene.lists.rds") 
  com.de.genelists=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/com.de.genelists.rds")
  
  GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
             "4k","5k","GSM4712885","GSM4712907")
  
  
  res.lists=list()
  for( i in 1 :1){
    GSM=GSMlists[[i]]
    
    if(i %in% c(1:4)){
      ref.gene<- ref.genelists[["Mouse"]]
      db=org.Mm.eg.db
    }else if( i %in% c(5:8) ){
      ref.gene<- ref.genelists[["Human"]]
      db= org.Hs.eg.db
    }
    
   
   com.gene=com.de.genelists[[GSM]]
    
   
   com.gene$gene=sub("_.*", "", com.gene$gene)
   
    loc=match(com.gene$gene,ref.gene$gene_id)
    com.gene$gene_name=ref.gene$gene_name[loc]
    
    if(anyNA(loc)==TRUE){
      stop(paste0(GSM," (",tool,")",":基因ensemble与symbol信息不完全匹配"))
    }
    
    #Gene enrichment analysis
    res <-enrichGO(gene =com.gene$gene_name,  
                        OrgDb = db, 
                        keyType ="SYMBOL",  
                        ont="ALL",         
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff =0.05,
                      
    ) 
    
    res.lists[[GSM]][["com.gene"]]=res@result
    
    
    
    tool.lists=names(deapa.gene.lists[[GSM]])
    
    for(j in 1: length(tool.lists)){
  
      tool=tool.lists[j]
      
      de.gene=deapa.gene.lists[[GSM]][[tool]]
      
      if(nrow(de.gene)==0){
        next
      }else{
        
        de.gene$gene_id=sub("_.*", "", de.gene$gene_id)
        
        loc=match(de.gene$gene_id,ref.gene$gene_id)
        de.gene$gene_name=ref.gene$gene_name[loc]
        
        if(anyNA(loc)==TRUE){
          stop(paste0(GSM," (",tool,")",":基因ensemble与symbol信息不完全匹配"))
        }
        
       
        temp.res <-enrichGO(gene =de.gene$gene_name,   
                            OrgDb = db, 
                            keyType ="SYMBOL",   
                            ont="ALL",        
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff =0.05,
                            
        ) 
        
        res.lists[[GSM]][[tool]]=temp.res@result
      
      }
      
      message(paste0(GSM,"(",tool,")",": ",j,"/",length(tool.lists)))
      
    }
    
  }
  
  
 #saveRDS(res.lists,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/old/go_res.lists.rds") 
  saveRDS(res.lists,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/go_res.lists.rds") 
  
 #(5.1)Bubble plot--------------------------------------------------------
 res.lists=readRDS("G:/scAPA/summary/temp/fig3/dexseq/go_res.lists.rds")
 
 GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
            "4k","5k","GSM4712885","GSM4712907")
 
 labs=c("Mouse_Sperm","Mouse_T_Cell","Mouse_TIP","Mouse_Stem_Cell",
  "Human_Pbmc2","Human_Pbmc3","Human_Covid_Pbmc1","Human_Covid_Pbmc2")
 
 
 GSM=GSMlists[1]
 lab=labs[1]
 #top5 terms
 tool.lists=names(res.lists[[GSM]])
 
 res=data.frame()
 for (j in 1: length(tool.lists)){
   tool=tool.lists[j]
   temp.res=subset(res.lists[[GSM]][[tool]],res.lists[[GSM]][[tool]]$ONTOLOGY=="BP")
  
   temp.res=temp.res[,c(1,3,4,5,7,10)]
   
   temp.res=temp.res %>%
     mutate(log10padj=-log10(p.adjust),
            method=tool)%>%
     arrange(p.adjust) %>%          
     slice_head(n = 5) 
   
   res=rbind(res,temp.res)
  
 }
 
 
 input.data=res
 
term.order <- res %>%
  group_by(Description) %>%
  summarise(
    count = n(),              
    min_padj = min(p.adjust)       
  ) %>%
  arrange(
    desc(count),             
    min_padj                  
  )
 
#tool order based on recall
data.res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.res.csv")

data.res$precision[is.na(data.res$precision)] <- 0


tool.order =data.res[data.res$sample==GSMlists[1],] %>% 
  group_by(tool)%>%
  mutate(mean.precision=mean(precision, na.rm = TRUE),
         mean.recall=mean(recall, na.rm = TRUE))%>%
  arrange(desc(mean.recall))


 input.data$Description<- factor(input.data$Description,
                              levels =term.order$Description)
 input.data$method <- factor(input.data$method,
                             levels =c("com.gene",tool.order$tool),
                             labels = c("Consensus Gene",tool.order$tool))
 
 
 
 
 #col <- colorRampPalette(c("#0f86a9","white","#fc8452"), alpha = TRUE)
 col <- colorRampPalette(c("#5A827E","#FCF8E8","#ECB390"), alpha = TRUE)
 
 p= ggplot(input.data, aes(x = method, y = Description, fill =log10padj)) + #fill =-1*log10(p.adj)
   geom_tile(aes(fill = rectheat),color="black",fill="white") + 
   geom_point(
     aes(size = Count,         
         fill =log10padj) ,        
         shape = 21 , 
     color = "#909fac",     
     stroke = 0.8          
   ) +
   
   scale_size_continuous(
     range = c(5, 10) ,    
     guide = guide_legend(         
       override.aes = list(shape = 21))
     #,guide = "none"         
   ) + 
   # 形状映射设置（显著=实心圆，不显著=星号）
   #scale_shape_manual(
   #  values = c("Not significant" = 8, "Significant" = 22),
   #  guide = guide_legend(title = "")
   #) +
   scale_fill_gradientn(
     colors =col(10))+  
   theme_minimal()+
   theme_bw()+
   labs(x=element_blank(),y=element_blank(),
        title =lab,
        fill="-log10(p.adjust)",
        size="Gene Number") + 
   theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
         axis.text.y = element_text(angle = 45, hjust = 1,size = 7),
         panel.background = element_rect(fill =  "gray100", color = 'white'),
         panel.grid.major = element_line(color = 'white', linetype = 'dotted'),
         panel.grid.minor = element_line(color = 'green', size = 2),
         legend.title = element_text(size = 7),  
         panel.grid = element_blank(),        
         axis.title.y = element_blank()
   ) 
 
  
 topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/go.res(",lab,").pptx"), width = 7, height = 7.5)
 
 #保存为png
 print(p)
 graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/go.res(",lab,").png") ,width=7,height=7.5)  
  
  
  
 ###5.Supplementary Figure----------------------------------------------------------------
 
##(1)Number of consensus DE APA genes for each sample-----------------------
 sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)

 res=read.csv("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.res.csv")
 loc=match(res$sample,sample.data$Sample_ID)
 res$species=sample.data$Species[loc]
 
 
 data=res[,c(8,1,3)]
 data=unique(data)
 
 data$species<-factor(data$species,
                      levels = c("Mouse","Human"))
 
 order=data %>%
   group_by(species)%>%
   arrange(desc(com.gene.num),.by_group = TRUE)
   
 
 
 data$sample <- factor(data$sample,
           levels=order$sample,
           labels = c("Mouse_Sperm","Mouse_Stem_Cell","Mouse_TIP","Mouse_T_Cell",
                      "Human_Covid_Pbmc2", "Human_Pbmc3","Human_Covid_Pbmc1","Human_Pbmc2"))
 
 col=c("Mouse" ="#FD8D3C","Human"='#163859',"Arabidopsis"="#064D4B")
 
 #col=c("Mouse" ="#fed6a2","Human"="#0980a3","Arabidopsis"="#fce9df")
p=ggplot(data, 
         aes(x=sample, y=com.gene.num,
             fill=species)) +
  geom_bar(
    stat = "identity"
    ,alpha=0.7
    #width = 0.7,                  
    #position = position_dodge(0.8) 
  )+
  
  geom_text(
    aes(
      label =com.gene.num ),  
    position = position_dodge(0.8),  
    vjust = -0.5,       
    size = 2.5            
  ) +
  p_theme+
  labs(x="",y=bquote("No. of Consensus Gene "),
       #fill="Species",
       title = "")+
  guides(fill = guide_legend(title = 'Species'))+
  scale_fill_manual(values=col)+
  theme(
    axis.text.x = element_text(angle = 45,hjust=1,size = 7)
  )


topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/con.gene.num.pptx"), width = 4, height = 3.5)

#保存为png
print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/con.gene.num.png") ,width=4,height=3.8) 
 
##(2)Agreement score for the DEAPA gene------------------------------------- 
deapa.gene.lists=readRDS("G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/deapa.gene.lists.rds")

sample.data <- read.csv("G:/scAPA/summary/temp/fig1/apa.number.csv",header = T)
GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848",
           "4k","5k","GSM4712885","GSM4712907")


 
res=data.frame()
for(i in 1:length(GSMlists)){
  GSM=GSMlists[i]
  species <- unique(sample.data[sample.data$Sample_ID==GSM,]$Species)
  protocol <- unique(sample.data[sample.data$Sample_ID==GSM,]$Protocol)
  
  temp=as.data.frame(table(unlist(lapply(deapa.gene.lists[[GSM]], function(df) df[, 1]))))
  
  temp$pa.num=nrow(temp)
  temp$tools.num=length(names(deapa.gene.lists[[GSM]]))
  #Agreement score
  temp$tools.score=(temp$Freq-1)/(temp$tools.num-1)
  temp$score=sum(temp$tools.score/temp$pa.num)
  temp$sample=GSM
  temp$protocol=protocol
  temp$species=species
  
  res=rbind(res,temp)
  
  }


write.csv(res,file="G:/scAPA/summary/temp/fig3/dexseq/based_com.gene/new/de.gene.agreement.score.csv",row.names = FALSE)  



###plot
res$species<-factor(res$species,
                     levels = c("Mouse","Human"))

data=unique(res[,c(6:9)]) %>%
     group_by(species) %>%
    arrange(desc(score),.by_group = TRUE)
      
data$sample<- factor(data$sample,
                     levels =data$sample,
                     labels = c("Mouse_Sperm","Mouse_TIP","Mouse_Stem_Cell","Mouse_T_Cell",
                                "Human_Covid_Pbmc2", "Human_Pbmc3","Human_Pbmc2","Human_Covid_Pbmc1"))




col=c("Mouse" ="#FD8D3C","Human"='#163859',"Arabidopsis"="#064D4B")

#col=c("Mouse" ="#fed6a2","Human"="#0980a3","Arabidopsis"="#fce9df")
p=ggplot(data, 
         aes(x=sample, y=score,
             fill=species)) +
  geom_bar(
    stat = "identity"
    ,alpha=0.7
    #width = 0.7,                 
    #position = position_dodge(0.8) 
  )+
  geom_text(
    aes(#label = sprintf("%.3f", score)),  
      label =round(score,3) ), 
    position = position_dodge(0.8),  
    vjust = -0.5,       
    size = 2.5           
  ) +
  p_theme+
  labs(x="",y=bquote("Agreement Score"),
       #fill="Species",
       title = "")+
  #scale_x_reverse() +  
  #scale_y_discrete(position = "right") +  
  #guides(fill='none')+  #隐藏图例
  guides(fill = guide_legend(title = 'Species'))+
  scale_fill_manual(values=col)+
  theme(
    axis.text.x = element_text(angle = 45,hjust=1,size = 7)
    
  )

topptx(p, filename = paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/de.gene.agreement.score.pptx"), width = 4, height = 3.5)

#保存为png
print(p)
graph2png(file=paste0("G:/scAPA/summary/Figure/fig3/dexseq/based_com.gene/new/de.gene.agreement.score.png") ,width=4,height=3.8) 

