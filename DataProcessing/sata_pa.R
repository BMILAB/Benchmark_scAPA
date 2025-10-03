
#===============================================================================
# 
#Goal:  Statistics predicted pA from different tools
#       (PACs filtered by lees expressed in 5% cells)
#
#
#                           Author: BXY
#                           Date:   20240830
#                           updata: 20250318
#===============================================================================
rm(list=ls())
options(stringsAsFactors = FALSE)
#library(liftOver)
library(rtracklayer)
library(plyr)
library(ggplot2)
library(extrafont)
library(export)
library(rtracklayer)
library(GenomicFeatures)
library(movAPA)
library(wesanderson)
library(dplyr)
library(patchwork)
library(aplot)
library(Rsamtools)

out.path <- "G:/scAPA/summary"
if(!dir.exists(out.path)){
  dir.create(out.path)
}
setwd(out.path)

###1.load data ( have runned all sample  data)------------------------------------------------------------------
###-------------------------------------------------
## PA table
# Each element of list is PA table(列表中每一个元素为一份数据，每份数据下又包含每个工具识别的pA结果)


##将同一样本数据集下不同工具识别的pA位点以及量化的pA表达量储存为列表并统计不同工具在基因组不同区域识别的pA位点数量

GSElists=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
           "arab")

all_datalists=c("merge","merge","merge","merge","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399","GSM2333586","GSM2333587",
            "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
            "GSM3490689","GSM3490690","GSM3490691")

##同一样本下不同数据集的样本编号（自定义）
all_GSMlists=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399",
               "GSM2333586","GSM2333587", "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
               "GSM3490689","GSM3490690","GSM3490691")

#protocols=c("10X","10X","10X","10X","10X","Microwell-seq","Drop-seq","CEL-seq",
#            "10X","10X","10X")

specieslists=c("Mouse","Human","Arabidopsis")

toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps","scutrquant")
tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps","scUTRquant")


#------------------------------------------------------------------------------------
#toolfile=toolfile.list[1]
#tool=tool.list[1]
#for(s in 1: length(GSElists)){
for(s in 9: 11){
  
sample=GSElists[s]
GSE=paste0("G:/scAPA/PAC/",sample)

species <- ifelse(s %in% 1:8, "Mouse",
             ifelse(s %in% 9:10, "Human",
                    ifelse(s == 11, "Arabidopsis",NA)))

protocol<- ifelse(s %in% c(1:5,9,10,11), "10x",
                  ifelse(s ==6, "Microwell-seq",
                         ifelse(s == 7, "Drop-seq",
                                ifelse(s == 8, "CEL-seq", NA))))

#file.path=paste0(GSE,"/",toolfile)
if(s %in% 1:4){
  ##同一样本下不同数据集的PAC数据存放文件编号
  datalists=c("merge")
  GSMlists=all_GSMlists[s]
  ##四份合并的数据集（mouse_sperm;mouse_tcell;mouse_tip:mouse_esc)的GSM编号均为merge
  #GSMlists=dir(file.path)
}else{
  GSMlists=dir(paste0(GSE,"/sierra"))
  datalists=GSMlists
}

if(sample=="mouse_bone"){
  toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps")
  tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps")
} else if(sample=="mouse_int"){
  toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","scraps")
  tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","scraps")
}else if(sample=="arab"){
  toolfile.lists=c("scapatrap","sierra","scape","scapa","maaper","polyapipe","scraps")
  tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","MAAPER","polyApipe","scraps")
}else{
  toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps","scutrquant")
  tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps","scUTRquant")
}


paclists=list()
utr.paclists=list()
sample.temp <- data.frame()

for(i in 1: length(GSMlists)){
  GSM=GSMlists[i]
  data=datalists[i]
  
  for(j in 1 :length(toolfile.lists)){
    message("start")
    
    toolfile=toolfile.lists[j]
    tool=tool.lists[j]
    message(paste0("read PACs from ",tool))
    
    ##未运行成功的数据
    conditions <- list(
      c(GSM == "GSM3490689", tool == "SCAPE"),
      c(GSM == "GSM2906399", tool == "scAPA"),
      c(GSM == "3k", tool == "scUTRquant"),
      c(GSM == "6k", tool == "scUTRquant")
    )
    
    if(any(sapply(conditions, all))){
      print("PAC data not exist")
    }else{
     
     
      if(tool == "scAPAtrap"){
        file.path=paste0(GSE,"/",toolfile,"/",data,"/APA.notails/PACds")
      }
      
      #no tail pac data,只有原始未过滤的数据，如需统计，去IP,并按表达量过滤后再重新注释
      if(tool == "scAPAtrap_1"){
        file.path=paste0(GSE,"/",toolfile,"/",data,"/APA.notails/PACds")
      }
      
      # tail pac data，但未除IP
      if(tool == "scAPAtrap_2"){
        file.path=paste0(GSE,"/",toolfile,"/",data,"/APA.tails/PACds")
      }
      
      if(!tool %in% c("scAPAtrap","scAPAtrap_1","scAPAtrap_2") ){
        file.path=paste0(GSE,"/",toolfile,"/",data,"/PACds")
      }
      scPACds=readRDS(paste0(file.path,"/anno.scPACds.rds"))
      paclists[[GSM]][[tool]]=scPACds
      
      scPACds3utr=readRDS(paste0(file.path,"/scPACds3utr.rds"))
      utr.paclists[[GSM]][[tool]]=scPACds3utr
      
      #创建数据框，记录每种工具在基因组不同区域识别的pa数量
      table(scPACds@anno$ftr)
      utr3=subset(scPACds@anno,scPACds@anno$ftr=="3UTR")
      utr5=subset(scPACds@anno,scPACds@anno$ftr=="5UTR")
      #cds=subset(scPACds@anno,scPACds@anno$ftr=="CDS")
      exon=subset(scPACds@anno,scPACds@anno$ftr=="exon")
      intron=subset(scPACds@anno,scPACds@anno$ftr=="intron")
      inter=subset(scPACds@anno,scPACds@anno$ftr=="intergenic")
      
      temp.data=data.frame(Species=species,
                           Sample=sample,
                           Sample_ID=GSM,
                           Protocol=protocol,
                           Tool=tool,
                           Total=nrow(scPACds@anno),
                           UTR3=nrow(utr3),
                           UTR5=nrow(utr5),
                           #CDS=nrow(cds),
                           Exon=nrow(exon),
                           Intron=nrow(intron),
                           Intergenic=nrow(inter)
                           
      )
      
      #temp.data=scPACds@anno %>%
      #  dplyr::group_by(ftr) %>%
      #  dplyr::summarise(count=n())
      #temp.data=as.data.frame(temp.data)
      #row.names(temp.data)=temp.data$ftr
      #temp.data=t(temp.data)
      #temp.data=temp.data[-1,]
      
      
      sample.temp <- rbind(sample.temp,temp.data)
      
    }
    
    
    message(paste0(sample,"(",GSM,")",": ",j,"/",length(toolfile.lists)))
        }
}
write.csv(sample.temp,file = paste0("G:/scAPA/summary/PACds/",sample,"_apa.number.csv"),row.names = FALSE)
saveRDS(paclists,file = paste0("G:/scAPA/summary/PACds/",sample,"_paclists.rds"))
saveRDS(utr.paclists,file = paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
#rm(paclists)
}


#-----------------------------------------------------------------------------------------------



##2.Convert to Granges object--------------------
toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps","scutrquant")
tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps","scUTRquant")


for(s in 1: 11){
  
  sample=GSElists[s]
  GSE=paste0("G:/scAPA/PAC/",sample)
  
  species <- ifelse(s %in% 1:8, "Mouse",
                    ifelse(s %in% 9:10, "Human",
                           ifelse(s == 11, "Arabidopsis",NA)))
  
  protocol<- ifelse(s %in% c(1:5,9,10,11), "10x",
                    ifelse(s ==6, "Microwell-seq",
                           ifelse(s == 7, "Drop-seq",
                                  ifelse(s == 8, "CEL-seq", NA))))
  #不同数据下不同工具的输出文件路径
  #file.path=paste0(GSE,"/",toolfile)
  if(s %in% 1:4){
    ##同一样本下不同数据集的PAC数据存放文件编号
    datalists=c("merge")
    GSMlists=all_GSMlists[s]
    ##四份合并的数据集（mouse_sperm;mouse_tcell;mouse_tip:mouse_esc)的GSM编号均为merge
    #GSMlists=dir(file.path)
  }else{
    GSMlists=dir(paste0(GSE,"/sierra"))
    datalists=GSMlists
  }
  
  if(sample=="mouse_bone"){
    toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps")
    tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps")
  } else if(sample=="mouse_int"){
    toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","scraps")
    tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","scraps")
  }else if(sample=="arab"){
    toolfile.lists=c("scapatrap","sierra","scape","scapa","maaper","polyapipe","scraps")
    tool.lists=c("scAPAtrap","Sierra","SCAPE","scAPA","MAAPER","polyApipe","scraps")
  }else{
    toolfile.lists=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps","scutrquant")
    tool.lists=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps","scUTRquant")
  }
  

gr_paclists=list()
utr3_gr_paclists=list()
for(i in 1: length(GSMlists)){
  GSM=GSMlists[i]
  data=datalists[i]
  
  for(j in 1 :length(toolfile.lists)){
    message("start")
    
    toolfile=toolfile.lists[j]
    tool=tool.lists[j]
    message(paste0("convert PACs from ",tool," to Granges"))
    
    
    ##未运行成功的数据
    conditions <- list(
      c(GSM == "GSM3490689", tool == "SCAPE"),
      c(GSM == "GSM2906399", tool == "scAPA"),
      c(GSM == "3k", tool == "scUTRquant"),
      c(GSM == "6k", tool == "scUTRquant")
    )
    
    if(any(sapply(conditions, all))){
      print("PAC data not exist")
    }else {
      if (tool =="scAPAtrap"){
      file.path=paste0(GSE,"/",toolfile,"/",data,"/APA.notails/PACds")
    } 
      if (tool == "scAPAtrap_1"){
      file.path=paste0(GSE,"/",toolfile,"/",data,"/APA.notails/PACds")
    } 
      if(tool == "scAPAtrap_2"){
      file.path=paste0(GSE,"/",toolfile,"/",data,"/APA.tails/PACds")
    }
      if (!tool %in% c("scAPAtrap","scAPAtrap_1","scAPAtrap_2") ){
        file.path=paste0(GSE,"/",toolfile,"/",data,"/PACds")}
      
    ###(1)ALL PAC----------------------------------------------------
    scPACds=readRDS(paste0(file.path,"/anno.scPACds.rds"))
    pac=scPACds@anno
    
    ##转换为grange对象
   
    #if(tool=="MAAPER" |tool=="SCAPE" | tool=="scraps" | tool=="scUTRquant"){
    #  pac=pac[,c( "chr", "strand", "coord","gene")]
    #  pac$start=pac$coord
    #  pac$end=pac$coord
     # pac=pac[,c( "chr", "strand", "start","end", "coord","gene")]
    #}else{
    #  pac=pac[,c( "chr", "strand", "start","end", "coord","gene")]
    #}
    
      pac=pac[,c( "chr", "strand", "coord","gene")]
      pac$start=pac$coord
      pac$end=pac$coord 
      
    
    #colnames(pa)=c( "chr", "strand", "start","end", "coord","gene")
    pac$count=Matrix::rowSums(scPACds@counts)
    pac=makeGRangesFromDataFrame(pac,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field=c("chr"),
                                 start.field="start",
                                 end.field="end",
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)
    
    gr_paclists[[GSM]][[tool]]=pac
    
    ##(2)3UTR PAC----------------------------------
    
    
    scPACds3utr=readRDS(paste0(file.path,"/scPACds3utr.rds"))
    utr_pac=scPACds3utr@anno
    
    ##转换为grange对象
    #if(tool=="MAAPER" | tool=="SCAPE" | tool=="scraps" | tool=="scUTRquant"){
    #  utr_pac=utr_pac[,c( "chr", "strand", "coord","gene")]
    #  utr_pac$start=utr_pac$coord
    #  utr_pac$end=utr_pac$coord
    #  utr_pac=utr_pac[,c( "chr", "strand", "start","end", "coord","gene")]
    #}else{
    #  utr_pac=utr_pac[,c( "chr", "strand", "start","end", "coord","gene")]
    #  }
      
    utr_pac=utr_pac[,c( "chr", "strand", "coord","gene")]
      utr_pac$start=utr_pac$coord
      utr_pac$end=utr_pac$coord
   
    #colnames(pa)=c( "chr", "strand", "start","end", "coord","gene")
    utr_pac$count=Matrix::rowSums(scPACds3utr@counts)
    utr_pac=makeGRangesFromDataFrame(utr_pac,
                                      keep.extra.columns=TRUE,
                                      ignore.strand=FALSE,
                                      seqinfo=NULL,
                                      seqnames.field=c("chr"),
                                      start.field="start",
                                      end.field="end",
                                      strand.field="strand",
                                      starts.in.df.are.0based=FALSE)
    
          utr3_gr_paclists[[GSM]][[tool]]=utr_pac
    }
          message(paste0(sample,"(",GSM,")",": ",j,"/",length(toolfile.lists)))
  }
}

saveRDS(gr_paclists,file = paste0("G:/scAPA/summary/PACds/",sample,"_gr_paclists(coord).rds"))
saveRDS(utr3_gr_paclists,file = paste0("G:/scAPA/summary/PACds/",sample,"_utr3_gr_paclists(coord).rds"))

}

###(3)mm10 refpa---------------------------------------------------------
##将真实pA位点的位置信息转换为Grange对象
#mm10的refPA
mouse_refpa=readRDS("G:/ref_file/refPAC/mouse/mm10_refPA.rds")


mouse_refpa=mouse_refpa[,c( "chr", "strand", "start","end", "coord","polyA_ID")]
#colnames(known.pa)=c( "chr", "strand", "start","end", "coord"," gene", "PAC_id")
#删除任意带NA值的行
#known.pa=na.omit(known.pa)

mouse_refpa=makeGRangesFromDataFrame(mouse_refpa,
                                  keep.extra.columns=TRUE,
                                  ignore.strand=FALSE,
                                  seqinfo=NULL,
                                  seqnames.field=c("chr"),
                                  start.field="start",
                                  end.field="end",
                                  strand.field="strand",
                                  starts.in.df.are.0based=FALSE)
#将mm39的真实pA替换为mm10的
refpa=readRDS("G:/ref_file/refPAC/collection_priorPACs_20240527.Rdata")
refpa[["Mouse"]]=mouse_refpa
saveRDS(refpa,file="G:/scAPA/summary/temp/fig1/new_refpa.rds")


##3.构建human/mouse/arabidopsis参考基因组序列列表-------------------------------
#(1)load genome--------------------------------------
fa.lists=list()

fafile<-"G:/ref_file/GRcm38fa+gtf/genecode/GRCm38.p6.genome.fa"
fa=FaFile(fafile)
indexFa(fa$path)
fa.lists[["Mouse"]]=fa

fafile<-"G:/ref_file/GRch38fa+gtf/v34/GRCh38.p13.genome.fa"
fa=FaFile(fafile)
indexFa(fa$path)
fa.lists[["Human"]]=fa

fafile<-"G:/ref_file/TAIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
fa=FaFile(fafile)
indexFa(fa$path)
fa.lists[["Arabidopsis"]]=fa

saveRDS(fa.lists,file = "G:/ref_file/fa.lists.rds")

##(2)load annotation genome---------------------------------------------
annotation.lists=list()

reference.file <- "G:/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf"
gff=makeTxDbFromGFF(reference.file , format="gtf")
gff<-parseGenomeAnnotation(gff)
annotation.lists[["Mouse"]]=gff

reference.file <- "G:/ref_file/GRch38fa+gtf/v34/gencode.v34.annotation.gtf"
gff=makeTxDbFromGFF(reference.file , format="gtf")
gff<-parseGenomeAnnotation(gff)
annotation.lists[["Human"]]=gff

reference.file <- "G:/ref_file/TAIR/Arabidopsis_thaliana.TAIR10.49.gtf"
gff=makeTxDbFromGFF(reference.file , format="gtf")
gff<-parseGenomeAnnotation(gff)
annotation.lists[["Arabidopsis"]]=gff

if(!dir.exists("G:/scAPA/summary/temp")){
  dir.create("G:/scAPA/summary/temp")
}

saveRDS(annotation.lists,file = "G:/scAPA/summary/temp/fig1/annotation.lists.rds")


###4.Annotated refpa------------------------------------------------------------
#source("C:/R data/果蝇节律神经元/APA差异分析/多源APA合并/APA combined/APA combined.function.R")

##加载refpa和参考基因组坐标
refPAC=readRDS("G:/scAPA/summary/temp/fig1/new_refpa.rds")
gff.lists<-readRDS("G:/scAPA/summary/temp/fig1/annotation.lists.rds")

species=c("Mouse","Human","Arabidopsis")
ext3utr=c(2000,2000,1000)
sources=c("GENCODE/PolyA_DB3","GENCODE/PolyA_DB3","ISO-seq/DRS")

anno.utr3.refpa=list()
anno.refpa=list()
anno.ref.pacds=list()
##注释refPAC的位点，获取PA位点所在基因
for(i in 1:length(species)){
  
  org=species[i]
  source=sources[i]
  cat(paste0("Annotate: ",org," refpA"))
  
  refpa=refPAC[[org]]
  refpa=as.data.frame(refpa)
  colnames(refpa)[1:3]=c("chr","UPA_start","UPA_end")
  refpa$coord=0
  #当strand="+",将识别到的peak的结束位置作为pA位点的坐标
  refpa[refpa$strand == "+",]$coord <- refpa[refpa$strand == "+",]$UPA_end
  #当srand="-",将识别到的peak的起始位置作为pA位点的坐标
  refpa[refpa$strand == "-",]$coord <- refpa[refpa$strand == "-",]$UPA_start
  refpa=refpa[,c("chr","UPA_start","UPA_end","strand","coord")]
  
  annoDB<-readPACds(refpa)
  annoDB<-annotatePAC(annoDB,gff.lists[[org]])
  annoDB<- ext3UTRPACds(annoDB,ext3UTRlen =ext3utr[i],extFtr='3UTR')
  anno.ref.pacds[[org]]=annoDB[annoDB@anno$ftr=="3UTR"]
  
  annoDB<-annoDB@anno
  annoDB<-subset(annoDB,!is.na(annoDB$gene))
  annoDB$PA_id <- paste0("PA:",annoDB$gene,":",annoDB$coord)

  annoDB$source=source
  
  annoDB <- annoDB[,c("PA_id","chr","UPA_start","UPA_end","strand","coord","ftr","gene","ftr_start","ftr_end","source")]
  #annoDB<-annoDB[,c("PA_id","chr","start","end","strand","coord","ftr","gene","ftr_start","ftr_end","source")]
  #annoMergePA需要gene_id这一列，所以需要将gene这一列改为相应的gene_id
  colnames(annoDB)<-c("PA_id","chr","start","end","strand","coord","ftr","gene_id","ftr_start","ftr_end","source")
  
  
  annoDB.3utr <- annoDB[annoDB$ftr == "3UTR",]
  
  
  head(annoDB)
  anno.refpa[[org]]=annoDB
  anno.utr3.refpa[[org]]=annoDB.3utr
  
}

saveRDS(anno.utr3.refpa,file = "G:/scAPA/summary/temp/fig1/anno.utr3.refpa.rds")
saveRDS(anno.refpa,file = "G:/scAPA/summary/temp/fig1/anno.refpa.rds")
saveRDS(anno.ref.pacds,file = "G:/scAPA/summary/temp/fig1/anno.utr3.ref.pacds.rds")

###4.1Select reference sequences for chi-square testing------------------------------------------------

##(1)计算refPA上下游100nt附近，碱基出现频率---------------------------------------------
bsgenome.list <- list()
bsgenome.list[["Human"]] <- BSgenome.Hsapiens.UCSC.hg38
bsgenome.list[["Mouse"]] <- BSgenome.Mmusculus.UCSC.mm10
bsgenome.list[["Arabidopsis"]] <- BSgenome.Athaliana.TAIR.TAIR9


refpa=readRDS("E:/scAPA/summary/temp/fig1/anno.utr3.ref.pacds.rds")
setwd("E:/scAPA/summary/refPA/fa")

species=c("Mouse","Human","Arabidopsis")

seq.num=1500
din.object.list <- list()
for( s in species){
  
  #org=species[s]
  ref.pac=refpa[[s]]
   
  if(s=="Arabidopsis"){
    ref.pac@anno$chr=paste0("Chr",ref.pac@anno$chr)
  }
  
  
  ##提取pA位点上下游100nt序列
  fafile = faFromPACds(ref.pac, bsgenome.list[[s]]
                        , what = "updn", fapre = paste0(s,"_","refPA.200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
  #读取DNA序列数据
  seq=readDNAStringSet(fafile, format="fasta")
  nseq=length(seq)
  if(nseq<seq.num){
    next
  }else{
    #计算序列集中每个位置上的核苷酸或氨基酸出现的频率，并生成一个矩阵，
    #其中列表示碱基，行表示序列集中的每个碱基在特定位置上出现的频率
    din<-consensusMatrix(seq, as.prob=TRUE, baseOnly=TRUE)
    din<- as.data.frame(t(din))
    din$other <- NULL
    din$pos <- c(1:nrow(din))
    din$number <- nseq
    #将宽数据转换为长数据格式，参数 id.vars，它表示在转换过程中要保持不变的列
    din.data2 <- reshape2::melt(din[,1:5],id.vars='pos',variable.name = 'base', value.name = 'freq')
    ATCGpicture.all <-ggplot(din.data2,aes(pos,y=freq,colour=base)) + geom_line() +
      xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw() +
      #scale_color_manual(values = base.col)+
      ggtitle("All-Whole-Gene-Body")
    
    din.object.list[[s]][["ref"]][["allseq"]] <- din
    din.object.list[[s]][["ref"]][["allfig"]] <- ATCGpicture.all
    
    #res=data.frame()
    for(i in 1:100){
      #循环100次，每循环一次从 nseq 个序列中随机抽取 5000/1500/1000 个不重复的序列。计算pA位点上下游100nt附近ATCG碱基出现的频率
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
        ggtitle(paste0(seq.num,"-Whole-Gene-Body"))
      
      din.object.list[[s]][["ref_sub"]][["seq"]][[i]] <- dinsub
      din.object.list[[s]][["ref_sub"]][["fig"]][[i]] <- ATCGpicture.sub
      
      
      #din.subdata$source.num=i
     #res=rbind(res,din.subdata)
      
    }
    
    ##计算100次结果中碱基平均分布频率
    #mean.res<- res %>%
    #  group_by(pos, base) %>%  # 按前两列分组
    #  summarise(
    #    mean.freq = mean(freq, na.rm = TRUE),  # 计算平均频率
    #    .groups = 'drop'  # 取消分组
    #  ) %>%
    #  arrange(pos)  # 按pos排序
    
    ###选取refPA随机的1000条序列的作为参考（这里第50次抽取的结果，没有特殊意义，就是单纯随机选了个结果）
    din.object.list[[s]][["ref"]][["repseq"]]=din.object.list[[s]][["ref_sub"]][["seq"]][[50]]
  }
}

print(din.object.list[["Mouse"]][["ref_sub"]][["fig"]][[50]])
din.object.list[["Mouse"]][["ref_sub"]][["seq"]][[50]][c(97:105),]

print(din.object.list[["Mouse"]][["ref"]][["allfig"]])
din.object.list[["Mouse"]][["ref"]][["allseq"]][c(97:105),]


#saveRDS(din.object.list,file="G:/scAPA/summary/refPA/refseq.1500.3UTR.rds")
saveRDS(din.object.list,file=paste0("G:/scAPA/summary/refPA/refseq.",seq.num,".3UTR.rds"))

##（2）计算所有参考序列与随机1500条序列间的卡方统计值，选值最小的序列作为最终的参考序列--------------------
din.object=readRDS("E:/scAPA/summary/refPA/refseq.1500.3UTR.rds")

species=c("Mouse","Human","Arabidopsis")

chiseq.object.list <-list()
res=data.frame()
for( s in species){
  
  #org=species[s]

  ###加载所有3UTR refPA的碱基分布频率
  ref.data <- din.object[[s]][["ref"]][["allseq"]]
  ref.data<- ref.data[,1:4]*ref.data$number
  ref.input <- c(ref.data[,1],
                 ref.data[,2],
                 ref.data[,3],
                 ref.data[,4])
  
  #重复抽样100次的结果
  temp=din.object[[s]][["ref_sub"]][["seq"]]
  for(n in 1:100){
    obj.data <- temp[[n]]
    obj.data<- obj.data[,1:4]*obj.data$number
    obj.input <- c(obj.data[,1],
                   obj.data[,2],
                   obj.data[,3],
                   obj.data[,4])
    #在1500个随机序列中的pA和所有refpA 上下游100nt附件ATCG出现的次数
    x <- matrix(c(as.integer(ref.input) ,
                  as.integer(obj.input)),nrow=2,byrow=TRUE)  
    
    chiseq.object.list[[s]][["statistic"]][[n]] <- chisq.test(x)$statistic
    chiseq.object.list[[s]][["pvalue"]][[n]] <- chisq.test(x)$p.value
    
    
  }
  
  ##整理100次抽样的结果
  
  statistic <- as.numeric(unlist(chiseq.object.list[[s]][["statistic"]]))
 
  pvalue <- as.numeric(unlist(chiseq.object.list[[s]][["pvalue"]]))
  
  temp.res <- data.frame(species=s,
                         statistic=statistic,
                         pvalue=pvalue)
  temp.res$id <- 1:100
  
  res=rbind(res,temp.res)
  
}

#statistic:201个碱基上的卡方统计值，最后区所有所有位置上的平均卡方统计值
res$label.value=res$statistic/201 

write.csv(res,file="E:/scAPA/summary/refPA/chiseq.sta.csv",row.names = FALSE)


#按species分组并找出每组label.value的最小值
data <- res %>%
  group_by(species) %>%
  filter(label.value == min(label.value)) %>%
  ungroup()

#选取卡方值最小的序列作为最终的参考序列
head(data)

din.object[["Mouse"]][["ref"]][["repseq"]]=din.object[["Mouse"]][["ref_sub"]][["seq"]][[87]]
din.object[["Human"]][["ref"]][["repseq"]]=din.object[["Human"]][["ref_sub"]][["seq"]][[16]]
din.object[["Arabidopsis"]][["ref"]][["repseq"]]=din.object[["Arabidopsis"]][["ref_sub"]][["seq"]][[79]]

saveRDS(din.object,file=paste0("E:/scAPA/summary/refPA/refseq.1500.3UTR.rds"))



##4.2Extract the 3UTR refpa and store it as a Granges object-------------------------------------
anno.utr3.refpa=readRDS("G:/scAPA/summary/temp/fig1/anno.utr3.refpa.rds")
anno.refpa=readRDS("G:/scAPA/summary/temp/fig1/anno.refpa.rds")

species=c("Mouse","Human","Arabidopsis")
utr3.refpa.list=list()
refpa.list=list()
for(i in 1:length(species)){
  
  org=species[i]
  cat(paste0("convert: ",org," refpA to Granges"))
  
  ###all pac
  refpa=anno.refpa[[org]]
  refpa=refpa[,c( "chr", "strand", "start","end", "coord","PA_id")]
  ###计算距离分布时，使用coord来计算，不使用peak region,故将其start和end替换为coord
  refpa$start=refpa$coord
  refpa$end=refpa$coord
  
      refpa=makeGRangesFromDataFrame(refpa,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("chr"),
                                     start.field="start",
                                     end.field="end",
                                     strand.field="strand",
                                     starts.in.df.are.0based=FALSE)
      refpa.list[[org]]=refpa
      
      ####3UTR pac
      
      utr3.refpa=anno.utr3.refpa[[org]]
      utr3.refpa=utr3.refpa[,c( "chr", "strand", "start","end", "coord","PA_id")]
      ###计算距离分布时，使用coord来计算，不使用peak region,故将其start和end替换为coord
      utr3.refpa$start=utr3.refpa$coord
      utr3.refpa$end=utr3.refpa$coord
      
      utr3.refpa=makeGRangesFromDataFrame(utr3.refpa,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=FALSE,
                                     seqinfo=NULL,
                                     seqnames.field=c("chr"),
                                     start.field="start",
                                     end.field="end",
                                     strand.field="strand",
                                     starts.in.df.are.0based=FALSE)
      
      
      utr3.refpa.list[[org]]=utr3.refpa
}
#saveRDS(utr3.refpa,file = "G:/scAPA/summary/temp/fig1/utr3.refpa.rds") 
saveRDS(utr3.refpa.list,file = "G:/scAPA/summary/temp/fig1/utr3.refpa(coord).rds")
saveRDS(refpa.list,file = "G:/scAPA/summary/temp/fig1/refpa(coord).rds")  

#5.Integrate data from all samples-----------------
###(1)single cell data------------------------------------------------------------

#Sample_ID=c("GSM2803334+GSM2803335","GSM2833284+GSM2833286","MI_DAY3_TIP+SHAM_TIP","GSM3629847+GSM3629848")
samples=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
          "arab")
#加载每个工具在基因组不同区域识别的pA数量统计表以及3UTR pA位点的grange对象
sample.data=data.frame()
utr.paclists=list()
gr_paclists=list()
utr3_gr_paclists=list()
for(i in 1: 11){
  sample=samples[i]
  
  data.path=paste0("G:/scAPA/summary/PACds/",sample,"_apa.number.csv")
  temp.data=read.csv(data.path,header = T)
  sample.data=rbind(sample.data,temp.data)
  
  ##scPACdataset,仅3UTR pA位点数据
  #temp=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr.paclists.rds"))
  #utr.paclists=c(utr.paclists,temp)
  
  ##Granges对象
  temp.paclists=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_gr_paclists(coord).rds"))
  gr_paclists=c(gr_paclists,temp.paclists)
  
  utr3.temp.paclists=readRDS(paste0("G:/scAPA/summary/PACds/",sample,"_utr3_gr_paclists(coord).rds"))
  utr3_gr_paclists=c(utr3_gr_paclists,utr3.temp.paclists)
  
  cat(paste0(sample,": ",i,"/",length(samples)))
}
#names(paclists)=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc_mef")# for four data
write.csv(sample.data,file = paste0("G:/scAPA/summary/temp/fig1/","apa.number.csv"),row.names = FALSE)
#saveRDS(utr.paclists,file = paste0("G:/scAPA/summary/temp/fig1/","utr.paclists.rds"))
saveRDS(gr_paclists,file = paste0("G:/scAPA/summary/temp/fig1/","gr_paclists(coord).rds"))
saveRDS(utr3_gr_paclists,file = paste0("G:/scAPA/summary/temp/fig1/","utr3_gr_paclists(coord).rds"))



###(2)bulk data---------------------------------------------------------------------
##pA位点信息来自原文或polyAsite数据库

samples=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc_mef")

paclists=list()
for(i in 1: length(samples)){
  sample=samples[i]
  
  ##PACdataset(only 3UTR)
  if(sample=="mouse_tip"){
    temp=readRDS(paste0("G:/scAPA/bulk/",sample,"/PACds/PACds3utr.rds")) 
  }else{
    temp=readRDS(paste0("G:/scAPA/bulk/",sample,"/RAW.PAC/PACds/PACds3utr.rds"))
  }
  
  paclists[[sample]]=temp
  
  
  cat(paste0(sample,": ",i,"/",length(samples)))
}

saveRDS(paclists,file = paste0("G:/scAPA/summary/temp/fig2/","bulk_paclists.rds"))


###来自polyAseqTrap
##mouse_sperm 样本
samples=c("2w","4w","5w")

paclists=list()
for(i in 1: length(samples)){
  sample=samples[i]
  path=paste0("G:/scAPA/bulk/mouse_sperm/polyaseqtrap/",sample,"/PACds")
  sralists=dir(path)
  for(j in 1:length(sralists)){
    ##PACdataset(only 3UTR)
    sra=sralists[j]
    file.path=paste0("G:/scAPA/bulk/mouse_sperm/polyaseqtrap/",sample,"/PACds/",sra,"/PACds3utr.rds")
    temp=readRDS(file.path)
    #添加样本标识符，方便整合
    colnames(temp@counts)=paste0(colnames(temp@counts),"_",sra)
    row.names(temp@colData)=colnames(temp@counts)
    paclists[[sra]]=temp
  }
  
  cat(paste0(sample,": ",j,"/",length(sralists)))
}

saveRDS(paclists,file = paste0("G:/scAPA/bulk/mouse_sperm/polyaseqtrap/merge/paclists.rds"))

##将不同样本识别的pA位点合并
pacds=mergePACds(paclists,d=24)

reference.file <- "G:/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf"#mouse
#gff<-parseGenomeAnnotation(reference.file)
gff=makeTxDbFromGFF(reference.file , format="gtf")
gff<-parseGenomeAnnotation(gff)
pacds<- annotatePAC(pacds,gff)

####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
pacds <- ext3UTRPACds(pacds,ext3UTRlen = 2000,extFtr='3UTR')
summary(pacds)
pacds<-pacds[pacds@anno$ftr=='3UTR']

saveRDS(pacds,file="G:/scAPA/bulk/mouse_sperm/polyaseqtrap/merge/PACds3utr.rds")

samples=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc_mef")


####********Create PAC lists for different sample*********#####################################
samples=c("mouse_sperm","mouse_tcell","mouse_tip")

paclists=list()
for(i in 1: length(samples)){
  sample=samples[i]
  
  ##PACdataset(only 3UTR)
  
  if (sample == "mouse_tip") {
    temp <- readRDS(paste0("G:/scAPA/bulk/", sample, "/PACds/PACds3utr.rds"))
  } else if (sample == "mouse_sperm") {
    temp <- readRDS(paste0("G:/scAPA/bulk/", sample, "/polyaseqtrap/merge/PACds3utr.rds"))
  } else if (sample == "mouse_tcell") {
    temp <- readRDS(paste0("G:/scAPA/bulk/", sample, "/RAW.PAC/PACds/PACds3utr.rds"))
  } else {
    stop("Unknown sample type!")  # 如果样本名不匹配，报错
  }

  
  paclists[[sample]]=temp
  
  cat(paste0(sample,": ",i,"/",length(samples)))
}

saveRDS(paclists,file = paste0("G:/scAPA/summary/temp/fig2/","bulk_paclists2.rds"))

