###read pA sites into PACdatase---------------------------------------------
##need R package
library(Sierra)
library(movAPA)
library(Seurat)
library(dplyr)
library(scAPA)
library(SCAPE)
library(dplyr)
library(polyApiper)
library(data.table)

##unify.pa---------------------------------------------------------------
##maaper
#anno="G:/scAPA/PAC/maaper/GSM2803334/pas.txt"
#pac=unify.pa(anno,tool="MAAPER")

##scapatrap
#pafile="G:/scAPA/PAC/scAPAtrap/GSM2803334/APA.tails/scAPAtrapData.rda"
#scPACds=unify.pa(pafile=pafile,tool="scAPAtrap")

unify.pa <- function(anno,count,pafile,tool){
  if(is.null(tool) | !tool %in% c("scAPAtrap","Sierra","SCAPTURE","scAPA","SCAPE","Infernape","MAAPER","polyApipe","scraps","scUTRquant","polyAseqTrap","QAPA","REPAC")){
    stop("tool should be scAPAtrap, Sierra, SCAPTURE, scAPA, SCAPE, Infernape,polyApipe, scraps MAAPER ,scUTRquant,polyAseqTrap,QAPA,REPAC ")
  }
  
  message(paste0("read PACs from ", tool))
  if(tool=="scAPAtrap"){
    PACds=read.scapatrap(pafile)
  }
  
  if(tool=="Sierra"){
    PACds=read.sierra(anno,count)
  }
  
  if(tool=="SCAPTURE"){
    PACds=read.scapture(anno,count) 
  }
  
  if(tool=="Infernape"){
   PACds=read.infernape(anno,count)
  }
  
  if(tool=="SCAPE"){
   PACds=read.scape(anno,count) 
  }
  
  if(tool=="scAPA"){
   PACds=read.scapa(anno,count)
  }
  
  if(tool=="MAAPER"){
  PACds=read.maaper(anno,count) 
  }
  
  
  if(tool=="polyApipe"){
    PACds=read.polyapipe(anno,count) 
  }
  
  if(tool=="scUTRquant"){
    PACds=read.scutrquant(pafile) 
  }
  
  if(tool=="scraps"){
    PACds=read.scraps(pafile) 
  }
  
  if(tool=="polyAseqTrap" ){
    PACds=read.bulk(pafile,method = tool)
  }
  
  if(tool=="QAPA" ){
    PACds=read.bulk(pafile,method = tool)
  }
  
  if(tool=="REPAC" ){
    PACds=read.bulk(pafile,method = tool)
  }
  
  return(PACds)
}







##scAPAtrap------------------------------------------------------------------

#pafile="G:/scAPA/PAC/scAPAtrap/GSM2803334/APA.tails/scAPAtrapData.rda"
#scPACds=read.scapatrap(pafile)
read.scapatrap <- function(pafile){
  load(pafile)
  scPACds=createPACdataset(counts=scAPAtrapData$peaks.count, anno=scAPAtrapData$peaks.meta)
  scPACds@colData$barcode<-row.names(scPACds@colData)
  #删掉group这一列
  #scPACds@colData<- subset(scPACds@colData, select= -group)
  scPACds@colData<- subset(scPACds@colData, select= barcode)
  
  print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells"))
  
  return(scPACds)
}

##sierra--------------------------------------------------------------------
#anno="G:/scAPA/PAC/sierra/GSM2803334/GSM2803334_peak_annotations.txt"
#count="G:/scAPA/PAC/sierra/GSM2803334/GSM2803334_peak.count"
#pac=read.sierra(anno,count)

read.sierra <- function(anno,count){
  
  peakfile <- anno
  countfile <- count
  
  peak_annotations <- read.delim(peakfile)
  #head(peak_annotations[1:3,])
  
  #Read in peak data saved in MEX format. Files can be in a gzipped (.gz) format.
  count <- ReadPeakCounts(data.dir=countfile,mm.file="matrix.mtx.gz",barcodes.file ="barcodes.tsv.gz")
  #library(Seurat)
  count <- as.data.frame(count)
  #head(count[1:3,1:4])
  
  #####保留count中与peak.annotion中相同的peak
  peak_annotations<-peak_annotations[intersect(row.names(count),row.names(peak_annotations)),]
  count<-count[intersect(row.names(count),row.names(peak_annotations)),]
  
  peak_annotations$coord <- 0
  
  peak_annotations[peak_annotations$strand == "+",]$coord <- peak_annotations[peak_annotations$strand == "+",]$end
  
  peak_annotations[peak_annotations$strand == "-",]$coord <- peak_annotations[peak_annotations$strand == "-",]$start
  
  anno <- peak_annotations[,c("seqnames","start","end","strand","coord","gene_id","pA_stretch","pT_stretch")]
  colnames(anno) <- c("chr","start","end","strand","coord","gene_name","pA_stretch","pT_stretch")
  
  anno <- unique(anno)
  
  ##筛选拥有gene name的PA位点
  anno<- anno[!is.na(anno$gene_name),]
  anno$gene_name <- gsub("[.].*","",anno$gene_name)
  count<-count[rownames(anno),]
  
  anno$polyA_ID<- rownames(anno)
  
  
  #combine the result of the annotated peak and peak count
  pafile<-cbind(anno,count)
  #head(pafile)
  coldata<-data.frame(barcode=colnames(pafile)[10:ncol(pafile)],row.names =colnames(pafile)[10:ncol(pafile)] )
  ##################读入PACdataset对象，进行pA位点的注释
  scPACds<-readPACds(pacFile= pafile,colDataFile =coldata,noIntergenic = FALSE,PAname = "PA")
  return(scPACds)
}

##scapture--------------------------------------------------------------------
#anno="G:/scAPA/PAC/mouse_sperm/scapture/GSM2803334_35/GSM2803334.PASquant.KeepPAS.metadata"
#count="G:/scAPA/PAC/mouse_sperm/scapture/GSM2803334_35/GSM2803334.PASquant.KeepCell.UMIs.tsv.gz"
#scPACds=read.scapture(anno,count)

read.scapture <- function(anno,count){
  #read PAS count matrix path
  count<- read.table(count,header=T)  
  
  rownames(count) <- gsub("_","-",gsub("\\|","-",count[,1],perl = TRUE),perl = TRUE)  
 
  colnames(count) <- sub(".1","-1",colnames(count))  
  #删除gene这一列
  count <- count[,-1]
  
  PASmetadata.raw<- data.frame(read.table(anno,header = FALSE))  #read PAS matadata file generated by SCAPTURE
  PASmetadata.raw$V4 <- gsub("_","-",gsub("\\|","-",PASmetadata.raw$V4,perl = TRUE),perl = TRUE) #rename PAS name from SCAPTURE 
  rownames(PASmetadata.raw) <- PASmetadata.raw$V4
  
  
  PASmetadata <- data.frame(polyA_ID = rownames(PASmetadata.raw))
  rownames(PASmetadata) <- PASmetadata$polyA_ID
  
  PASmetadata$gene_name <- as.character(PASmetadata.raw$V13)
  PASmetadata$chr <- as.character(PASmetadata.raw$V1)
  PASmetadata$start <- as.numeric(PASmetadata.raw$V2)
  PASmetadata$end <- as.numeric(PASmetadata.raw$V3)
  PASmetadata$strand <- as.character(PASmetadata.raw$V6)
  
  PASmetadata=unique(PASmetadata)
 
  PASmetadata<- PASmetadata[!is.na(PASmetadata$gene_name),]
  #保留两者共有的PAS
  PASmetadata<- PASmetadata[intersect( rownames(PASmetadata),rownames(count)),]
  count <- count[intersect( rownames(PASmetadata),rownames(count)),]
  
  #####添加PAS的坐标
  PASmetadata$coord <- 0                              
  PASmetadata[PASmetadata$strand == "+",]$coord <- PASmetadata[PASmetadata$strand == "+",]$end
  PASmetadata[PASmetadata$strand == "-",]$coord <- PASmetadata[PASmetadata$strand == "-",]$start
  
  scPACds <- createPACdataset(counts=count, anno=PASmetadata)
  scPACds@colData$barcode<-row.names(scPACds@colData)
  #删掉group这一列
  scPACds@colData<- subset(scPACds@colData, select= barcode)
  print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
  
  return(scPACds)
  
}


##Infernape--------------------------------------------------------------------
#anno="G:/scAPA/PAC/infernape/GSM2803334/anno_filtered.csv"
#count="G:/scAPA/PAC/infernape/GSM2803334/cnt_mat"
#scPACds=read.infernape(anno,count)

read.infernape <- function(anno,count=NULL){
  #read peak count
  if(!is.null(count)){
    counts.dir=count
    count = Matrix::readMM(paste0(counts.dir, '/matrix.mtx'))
    colnames(count) = unlist(utils::read.table(paste0(counts.dir, '/barcodes.tsv'), sep = '\t', header = F))
    rownames(count) = unlist(utils::read.table(paste0(counts.dir, '/sitenames.tsv'), sep = '\t', header = F))
    count <- as.data.frame(count)
    #head(cnt)
    ##read filtered peak annotation table
    anno.file=anno
    anno = utils::read.csv(anno.file, stringsAsFactors = FALSE)
    #head(anno)[,1:15]
    
    rownames(anno) = anno$polyA_ID
    #####保留count中与peak.annotion中相同的peak
    anno<-anno[intersect(row.names(count),row.names(anno)),]
    count<-count[intersect(row.names(count),row.names(anno)),]
    
    anno$coord <- 0
   
    anno[anno$strand == "+",]$coord <- anno[anno$strand == "+",]$to
   
    anno[anno$strand == "-",]$coord <- anno[anno$strand == "-",]$from
    
    #mode.pos:对原始计数应用高斯核平滑并生成降噪曲线。然后将峰值模态确定为该曲线的局部最大值
    anno <- anno[,c("seq","from","to","strand","coord","gene","mode.pos")]
    colnames(anno) <- c("chr","start","end","strand","coord","gene_name","mode.pos")
   
    anno <- unique(anno)
    
    ##筛选拥有gene name的PA位点
    anno<- anno[!is.na(anno$gene_name),]
    anno$gene_name <- gsub("[.].*","",anno$gene_name)
    count<-count[rownames(anno),]
    
    anno$polyA_ID<-rownames(anno)
    
    scPACds <- createPACdataset(counts=count, anno=anno)
    scPACds@colData$barcode<-row.names(scPACds@colData)
    #删掉group这一列
    scPACds@colData<- subset(scPACds@colData, select= barcode)
    print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
    
  }else{
    
    ##read filtered peak annotation table
    anno.file=anno
    anno = utils::read.csv(anno.file, stringsAsFactors = FALSE)
    #head(anno)[,1:15]
    
    rownames(anno) = anno$polyA_ID
    
    anno$coord <- 0
   
    anno[anno$strand == "+",]$coord <- anno[anno$strand == "+",]$to
    
    anno[anno$strand == "-",]$coord <- anno[anno$strand == "-",]$from
    
    #mode.pos:对原始计数应用高斯核平滑并生成降噪曲线。然后将峰值模态确定为该曲线的局部最大值
    anno <- anno[,c("seq","from","to","strand","coord","gene","mode.pos")]
    colnames(anno) <- c("chr","start","end","strand","coord","gene_name","mode.pos")
    #根据peak位置信息去除重复的peak数据
    anno <- unique(anno)
    
    ##筛选拥有gene name的PA位点
    anno<- anno[!is.na(anno$gene_name),]
    anno$gene_name <- gsub("[.].*","",anno$gene_name)
    anno$polyA_ID<-rownames(anno)
    
    scPACds <-readPACds(pacFile=anno)
    print(paste0(length(row.names(scPACds@anno))," PACs"))
  }
  
  
  return(scPACds)
}

##SCAPE-----------------------------------------------------------------------
#anno="G:/scAPA/PAC/arab/scape/GSM3490690/collapse_pa.tsv.gz"
#count="G:/scAPA/PAC/arab/scape/GSM3490690/pasite.csv.gz"
#scPACds=read.scape(anno,count)

read.scape <- function(anno,count){
  ##read peak count
 
  names(count) <- basename(dirname(count))
  
  ##read peak location infomation
  # load the collapse pA site file.
  # generate from `script/group_pa.py`
  collapse_pa <- anno
  
  #返回pa表达矩阵
  pa_count <- SCAPE::loadData(
    fileList = count,
    collapsePa = collapse_pa,
    matrix = TRUE,
    cores =1 
  )

  anno<- data.table::fread(collapse_pa)
  anno=as.data.frame(anno)
  label=unique(anno$label)
  colnames(pa_count)= gsub(paste0(label,"."),"",colnames(pa_count))
  #head(anno)
  anno=anno[,c("chrom","strand","collapse_pa")]
  colnames(anno)=c("chr","strand","coord")
  
  
  anno=unique(anno)
  anno$polyA_ID=paste0(anno$chr,":",anno$coord,":",anno$strand)
  row.names(anno)=anno$polyA_ID
  
  #####保留count中与annotion中相同的peak
  anno<-anno[intersect(row.names(pa_count),row.names(anno)),]
  count<-pa_count[intersect(row.names(pa_count),row.names(anno)),]
 
  scPACds <- createPACdataset(counts=count, anno=anno)
  scPACds@colData$barcode<-row.names(scPACds@colData)
  #删掉group这一列
  scPACds@colData<- subset(scPACds@colData, select= barcode)
  print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
  
  return(scPACds)
  
}

##scAPA------------------------------------------------------------------------
#count="G:/scAPA/PAC/mouse_bone/scAPA/GSM2906396/outs/cell_Peaks.RDS"
#anno="G:/scAPA/PAC/mouse_bone/scAPA/GSM2906396/outs/out_peaks.txt"
#scPACds=read.scapa(anno,count)

read.scapa <- function(anno,count){
  pa=readRDS(count)
  
  anno.count=pa@row.Data
  anno=read.table(anno, header = T, sep = "\t", stringsAsFactors = FALSE)
  

  dup=length(anno.count$GeneID)==length(anno$GeneID)
  if(dup==FALSE){
   
    anno.count=tidyr::separate(anno.count, Chr, into = c("Chr1", "Chr2"), sep = ";")
    anno.count=tidyr::separate(anno.count, Strand, into = c("Strand1", "Strand2"), sep = ";")
   
    anno.count=tidyr::separate(anno.count, Start, into = c("Start1", "Start2"), sep = ";")
    
    anno.count=tidyr::separate(anno.count, End, into = c("End1", "End2"), sep = ";")
    
   
    anno.count$chr=anno.count$Chr1
    anno.count$strand=anno.count$Strand1
    
    anno.count$start=pmin(anno.count$Start1,anno.count$Start2)
    anno.count[which(is.na(anno.count$start)),]$start=anno.count[which(is.na(anno.count$start)),]$Start1
    
    
    anno.count$end=pmax(anno.count$End1,anno.count$End2)
    anno.count[which(is.na(anno.count$end)),]$end=anno.count[which(is.na(anno.count$end)),]$End1
    
      row.names(anno.count)=anno.count$GeneID
      #which(duplicated(anno.count[,c("chr","strand","start","end")]))
      anno.count=anno.count[,c("chr","strand","start","end")]
      anno=anno.count
    
  }else{
    row.names(anno.count)=anno.count$GeneID
    anno.count=anno.count[,c("Chr","Strand","Start","End")]
    colnames(anno.count)=c("chr","strand","start","end")
    anno=anno.count
  }
  
  anno$coord <- 0
  
  anno[anno$strand == "+",]$coord <- anno[anno$strand == "+",]$end
  
  anno[anno$strand == "-",]$coord <- anno[anno$strand == "-",]$start
  #根据pa位置信息去除重复的pa数据
  anno=unique(anno)
  anno$polyA_ID=row.names(anno)
  
  ##exclude internal priming suspected peaks, peaks having a stretch of at least 8 consecutive As in the region between 10 nt to 140 nt. to the peak's end are filtered
  pa=filter_IP_cell(x = pa, int.priming.seq = "AAAAAAAA",
            left = 10,
            right = 140)
  #anno=anno[pa@row.Data$GeneID,]
  
  #@cells.counts储存pa在单个细胞中的表达量
  count=pa@cells.counts
  count=count[,-1]
  
  #####保留count中与annotion中相同的peak
  anno<-anno[intersect(row.names(count),row.names(anno)),]
  count<-count[intersect(row.names(count),row.names(anno)),]
  
  scPACds <- createPACdataset(counts=count, anno=anno)
  scPACds@colData$barcode<-row.names(scPACds@colData)
  #删掉group这一列
  scPACds@colData<- subset(scPACds@colData, select= barcode)
  print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
  
  return(scPACds)
    
}

##maaper----------------------------------------------------------------------
#anno="G:/scAPA/PAC/mouse_sperm/maaper/GSM2803334/pas.txt"
#count="G:/scAPA/PAC/mouse_sperm/maaper/GSM2803334/gene.txt"
#scPACds=read.maaper(anno,count)

read.maaper <- function(anno,count){
 
  count=read.table(count,header = T, sep = "\t", stringsAsFactors = FALSE,quote = "")
  
  anno=read.table(anno,header = T, sep = "\t", stringsAsFactors = FALSE,quote = "")
  
  loc=match(anno$gene,count$gene)
  anno$nread.c1=count$nread.c1.1[loc]
  anno$nread.c2=count$nread.c2.1[loc]
  
  anno$count.c1=round(anno$frac.c1*anno$nread.c1,0)
  anno$count.c2=round(anno$frac.c2*anno$nread.c2,0)
 
  
  #若count.c1==count.c2,则表示仅使用同一条件下的bam文件用于pA识别,count.c1=count.c2,保留一列即可
  if(identical( anno$count.c1, anno$count.c2)==TRUE){
    anno$count= anno$count.c1
    anno= subset(anno,select=-c(count.c1,count.c2))
    anno$polyA_ID=anno$pas
    anno=tidyr::separate(anno, pas, into = c("chr", "coord","strand"), sep = ":")
    anno=anno[,c("chr","coord","strand","gene","polyA_ID","count")]
    
    colnames(anno)=c("chr","coord","strand","gene_name","polyA_ID","count")
    
    #去除重复的peak数据
    anno <- anno[!duplicated(anno$polyA_ID),]
    
    ##筛选拥有gene name的PA位点
    anno<- anno[!is.na(anno$gene_name),]
    count=data.frame(row.names = anno$polyA_ID,count=anno$count)
    
    row.names(anno)=anno$polyA_ID
    anno=subset(anno,select = -count)
    
    scPACds<-createPACdataset(anno= anno,counts = count)
    
  } else{
    anno$polyA_ID=anno$pas
    anno=tidyr::separate(anno, pas, into = c("chr", "coord","strand"), sep = ":")
    anno=anno[,c("chr","coord","strand","gene","polyA_ID","count.c1","count.c2")]
    
    colnames(anno)=c("chr","coord","strand","gene_name","polyA_ID","count.c1","count.c2")
    
    #去除重复的peak数据
    anno <- anno[!duplicated(anno$polyA_ID),]
    
    ##筛选拥有gene name的PA位点
    anno<- anno[!is.na(anno$gene_name),]
    count=data.frame(row.names = anno$polyA_ID,count.c1=anno$count.c1,count.c2=anno$count.c2)
    
    row.names(anno)=anno$polyA_ID
    anno=subset(anno,select = -c(count.c1,count.c2))
    
    scPACds<-createPACdataset(anno= anno,counts = count)
  }
  
  
  print(paste0(length(row.names(scPACds@anno))," PACs"))
  
  return(scPACds)
  
}

###scutrquant（finish）-----------------------------------------------

#pafile=readRDS("G:/scAPA/PAC/mouse_esc/scutrquant/merge/merge.txs.Rds")

read.scutrquant <- function(pafile){
  pac=readRDS(pafile)
  #pac:SingleCellExperiment object. The counts in the object is a sparse Matrix of 3' UTR isoform counts; 
 
  ##pac@assays@data@listData$counts:转录本的表达信息
  count=pac@assays@data@listData$counts
 
  colnames(count)=gsub(".*_", "", colnames(count))
  
 
  raw.anno=as.data.frame(pac@rowRanges)
  raw.anno= raw.anno[,c("seqnames","start","end","width","strand","transcript_id")]
  colnames(raw.anno)[c(1,2,3,6)]=c("chr","truncation_start","truncation_end","group_id")
  
 
  anno<-  raw.anno %>%
    group_by(group_id) %>%
    mutate(coord = ifelse(strand == "+", max(truncation_end), min(truncation_start))) %>%
    filter((strand == "+" & truncation_end == coord) | (strand == "-" & truncation_start == coord)) %>%
    ungroup()
  
  
  anno$polyA_ID=paste0(anno$chr,":",anno$coord,":",anno$strand)
  #去除重复的peak数据
  anno <- anno[!duplicated(anno$polyA_ID),]
  anno=as.data.frame(anno)
  row.names(anno)=anno$group_id
  
  #保留相同的peak
  anno=anno[intersect(anno$group_id,row.names(count)),]
  count=count[intersect(anno$group_id,row.names(count)),]
  
  scPACds <- createPACdataset(counts=count, anno=anno)
  scPACds@colData$barcode<-row.names(scPACds@colData)
  #删掉group这一列
  scPACds@colData<- subset(scPACds@colData, select= barcode)
  print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
  
  return(scPACds)
}


##polyapipe---------------------------------------------------------------
#anno="G:/scAPA/PAC/mouse_tcell/polyapipe/merge/merge_polyA_peaks.gff"
#count="G:/scAPA/PAC/mouse_tcell/polyapipe/merge/merge_counts.tab.gz"
#scPACds=read.polyapipe(anno,count)

read.polyapipe <- function(anno, count){
  #读入表达矩阵
  the_counts <- read.table(count, sep="\t", header=TRUE, 
                           as.is=TRUE, 
                           colClasses = c("factor","factor", "numeric"),
                           nrows=-1)
  colnames(the_counts) <- c("peak","cell","count")
  count<- Matrix::sparseMatrix(
    i = as.integer(the_counts$peak),
    j = as.integer(the_counts$cell),
    x = as.integer(the_counts$count))
  colnames(count) <- levels(the_counts$cell)
  rownames(count) <- levels(the_counts$peak)   
  
  message("Loaded ", nrow(count), " x ", ncol(count), " matrix of counts")
  
  #读入pA位点位置信息
  anno=polyApiper::read_polyA_peak_file_gtf(anno)
  anno$coord <- 0
  
  anno[anno$strand == "+",]$coord <- anno[anno$strand == "+",]$end
  
  anno[anno$strand == "-",]$coord <- anno[anno$strand == "-",]$start
  
  anno$polyA_ID=paste0(anno$chr,":",anno$coord,":",anno$strand)
  #去除重复的peak数据
  anno <- anno[!duplicated(anno$polyA_ID),]
  anno =anno[,c("chr","start","end","strand","coord","peak","misprime")]
  row.names(anno)=anno$peak
  
  #保留相同的peak
  anno=anno[intersect(anno$peak,row.names(count)),]
  count=count[intersect(anno$peak,row.names(count)),]
  
  scPACds <- createPACdataset(counts=count, anno=anno)
  scPACds@colData$barcode<-row.names(scPACds@colData)
  #删掉group这一列
  scPACds@colData<- subset(scPACds@colData, select= barcode)
  print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
  
  return(scPACds)
  
}

###scraps----------------------------------------------------------
#pafile="G:/scAPA/PAC/mouse_gfp/scraps/MI_DAY7_GFP/MI_DAY7_GFP_R2_counts.tsv.gz"
#scPACds=read.scraps(pafile)


read.scraps <- function(pafile){
  #count table : +-10 around PolyA_DB sites, by cell barcode
 count <- scrapR::scraps_to_matrix(pafile, 
                                 n_min = 1,
                                 gene_min = 1,
                                 alt_only = FALSE,
                                 filter_low = FALSE,
                                 cell_ids = NULL,
                                 types = NULL,
                                 types2 = NULL,
                                 pf = NULL)
 
 anno=data.frame(row.names = row.names(count),peak=row.names(count))
 #anno <- tidyr::separate(anno, col = peak, into = c("gene_symbol", "number", "gene_name","chr","coord","strand","type"), sep = "[;_]")
 

 anno <- if (any(grepl(":", anno$peak))) {
   tidyr::separate(anno, col = peak, 
                   into = c("gene_name", "chr", "coord", "strand"), 
                   sep = ":")
 } else {
   tidyr::separate(anno, col = peak, 
                   into = c("gene_symbol", "number", "gene_name", "chr", "coord", "strand", "type"), 
                   sep = "[;_]")
 }
 
 
 anno$polyA_ID=paste0(anno$chr,":",anno$coord,":",anno$strand)
 #去除重复的peak数据
 anno <- anno[!duplicated(anno$polyA_ID),]
 
 #anno=anno[c("chr","coord","strand","gene_symbol","polyA_ID")]
 if ("gene_symbol" %in% colnames(anno)) {
   anno <- anno[c("chr", "coord", "strand", "gene_symbol", "polyA_ID")]
 } else {
   anno <- anno[c("chr", "coord", "strand", "gene_name", "polyA_ID")]
 }
 
 #保留相同的peak
 anno=anno[intersect(row.names(anno),row.names(count)),]
 count=count[intersect(row.names(anno),row.names(count)),]
 
 #id=match(row.names(count),row.names(anno))
 row.names(anno)=anno$polyA_ID
 row.names(count)=anno$polyA_ID
 
 
 scPACds <- createPACdataset(counts=count, anno=anno)
 scPACds@colData$barcode<-row.names(scPACds@colData)
 #删掉group这一列
 scPACds@colData<- subset(scPACds@colData, select= barcode)
 print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
 return(scPACds)
}






##3sep----------------------------------------------------
#pafile="G:/scAPA/bulk/mouse_esc/SRR1033836_PAtable_renew.Rdata"
#PACds=read.bulk(pafile,method="polyAseqTrap")

#pafile="G:/AML/REPAC/filter.se.rds"
#PACds=read.bulk(pafile,method="REPAC")

read.bulk <- function(pafile, method){
  if(method=="polyAseqTrap"){
   
    pac=readRDS(pafile)
    pac=pac[["pa.coord"]][["pa.coord"]]
    pac=pac[,c("PAid","seqnames","start","end","strand","coord","total.count")]
    
    anno=pac[,c("PAid","seqnames","start","end","strand","coord")]
    row.names(anno)=anno$PAid
    colnames(anno)[2]="chr"
    
    count=data.frame(row.names = pac$PAid,count=pac$total.count)
    
    PACds<-createPACdataset(anno= anno,counts = count)
    
  }
  
  if(method=="QAPA"){
    pac=read.delim(pafile)
    
    tpm.col <- grep("TPM", names(pac),value = TRUE)
    pac=pac[,c("Chr","UTR3.Start","UTR3.End","Gene_Name","APA_ID","Strand",tpm.col)]
    colnames(pac)[1:6]=c("chr","utr3.start","utr3.end","gene_name","paid","strand")
    
    pac$coord=0
    pac[pac$strand=="+",]$coord=pac[pac$strand=="+",]$utr3.end
    pac[pac$strand=="-",]$coord=pac[pac$strand=="-",]$utr3.start
    #过滤不含有“gene_name"的pA
    pac=pac[!is.na(pac$gene_name),]
    row.names(pac)=paste0(pac$chr,":",pac$gene_name,":",pac$coord,":",pac$strand)
    
    anno=pac[,c("chr","utr3.start","utr3.end","gene_name","paid","strand","coord")]
    
    if(length(tpm.col)==1){
      count=data.frame(row.names = row.names(pac),
                       tpm=pac[,tpm.col])
      colnames(count)=tpm.col
    }else{
      count=pac[,tpm.col]
      #colnames(count)=gsub("\\.","_",colnames(count))
    }
    
    PACds<-createPACdataset(anno= anno,counts = count)
    
    
    PACds=subsetPACds(PACds, totPACtag=1,minExprConds=0,verbose =T)
  }
  
  if(method=="REPAC"){
    pac=readRDS(pafile)
    anno=as.data.frame(pac@rowRanges)
    count=pac@assays@data$counts
    
    #过滤不含有“gene_name"的pA
    anno=anno[!is.na(anno$gene_name),]
    count=count[row.names(anno),]
    #colnames(anno)
    anno=anno[,c("seqnames","start","end","coord","strand","gene_name","geneID")]
    colnames(anno)=c("chr","start","end","coord","strand","gene_name","paid")
    
    PACds<-createPACdataset(anno= anno,counts = count)
    
    PACds=subsetPACds(PACds, totPACtag=3,minExprConds=2,verbose =T)
  }
  
  print(paste0(length(row.names(PACds@anno))," PACs"))
  return(PACds)
}








############################################################


annoMergePA <- function(pa_anno,annoDB,source_name,obs=10){
  #实现程序进度条，样式3使用竖线“|”指定进度区间，同时在进度条右侧显示百分进度
  pb <- txtProgressBar(style=3)
  star_time <- Sys.time()  
  
  pa_anno$PA_id <- paste0("PA:",pa_anno$gene,":",pa_anno$coord)
  pa_anno$PA_id2 <- pa_anno$PA_id
  pa_anno$new.coord<-pa_anno$coord
  #pa_anno$gene <- toupper(pa_anno$gene)#toupper:将小写字符串转换为大写字符串
  
  #annoDB$gene_id <-  toupper(annoDB$gene_id)
  sub_anno <- pa_anno[(pa_anno$gene %in% annoDB$gene_id),]
  
  for(i in 1:nrow(sub_anno)){
    setTxtProgressBar(pb, i/nrow(sub_anno))
    
    gene_id <- sub_anno[i,]$gene
    coord <- sub_anno[i,]$coord
    start <- sub_anno[i,]$start
    end <- sub_anno[i,]$end
    ftr <- sub_anno[i,]$ftr
    PA_id <- sub_anno[i,]$PA_id
    ftr_start <- sub_anno[i,]$ftr_start
    ftr_end <- sub_anno[i,]$ftr_end
    
    sub_db <- subset(annoDB,annoDB$gene_id == gene_id)
    #sub_db <- annoDB[annoDB$gene_id == gene_id,]
    ref_start <- sub_db$start
    ref_end <- sub_db$end
    
    
    #判断是否存在区间交集
    #obs=10表示一个范围误差
    #排序
    ranks <- rank(ref_start)
    ref_start <- ref_start[ranks] -obs
    ref_end <- ref_end[ranks] +obs
    #loc <- start < ref_end & end > ref_start
    loc <- ref_start < end & ref_end > start
    if(sum(loc) == 1){
      #当pa位点的起始和结束区间与参考坐标系的某些位点只存在一个重合区间时，即有唯一值
      sub_anno[i,]$PA_id2 <- sub_db[ranks,]$PA_id[loc]
      sub_anno[i,]$new.coord<-sub_db[ranks,]$coord[loc]
    }else if(sum(loc) > 1){
      #当有多个重合区间时，通过距离计算筛选出唯一值，用参考坐标系中距离最近的pA位点进行注释
      ref_coord <- as.numeric(sub_db[ranks,]$coord[loc])
      loc2 <- which.min(abs(ref_coord-coord))
      sub_anno[i,]$PA_id2 <- sub_db[ranks,]$PA_id[loc][loc2]
      sub_anno[i,]$new.coord <- sub_db[ranks,]$coord[loc][loc2]
      
    }else{
      dt <- sub_db[1,]
      dt$PA_id <- PA_id
      dt$start <- start
      dt$end <- end
      dt$coord <- coord
      dt$ftr_start <- ftr_start
      dt$ftr_end <- ftr_end
      dt$source <- source_name
      annoDB <- rbind(annoDB,dt)
    }
  }
  
  #返回原始的anno
  locs <- match(rownames(sub_anno),rownames(pa_anno))
  pa_anno$PA_id2[locs] <- sub_anno$PA_id2
  pa_anno$new.coord[locs] <- sub_anno$new.coord
  
  #DB整合进不包含的PA
  
  No_anno <- pa_anno[!(pa_anno$gene %in% annoDB$gene_id),]
  new_anno <- No_anno[,c("PA_id","chr","start","end","strand","coord","ftr","gene","ftr_start","ftr_end")]
  
  if(nrow(new_anno) >1){
    new_anno$source <- source_name
    colnames(new_anno)<- colnames(annoDB)
    annoDB <- rbind(annoDB,new_anno)
  }
  
  
  ##
  end_time <- Sys.time() 
  close(pb)
  run_time <- end_time - star_time
  
  return(list(pa_anno,annoDB))
  
}

#######################################################将各样本数据的pa位点merge后再根据参考坐标系进行重新注释
annoMergePA.2 <- function(pa_anno,annoDB,source_name,obs=10){
  #实现程序进度条，样式3使用竖线“|”指定进度区间，同时在进度条右侧显示百分进度
  pb <- txtProgressBar(style=3)
  star_time <- Sys.time()  
  
  pa_anno$PA_id <- paste0("PA:",pa_anno$gene,":",pa_anno$coord)
  pa_anno$PA_id2 <- pa_anno$PA_id
  
  pa_anno$gene <- toupper(pa_anno$gene)#toupper:将小写字符串转换为大写字符串
  
  annoDB$gene <-  toupper(annoDB$gene)
  sub_anno <- pa_anno[(pa_anno$gene %in% annoDB$gene),]
  
  for(i in 1:nrow(sub_anno)){
    setTxtProgressBar(pb, i/nrow(sub_anno))
    
    gene_id <- sub_anno[i,]$gene
    coord <- sub_anno[i,]$coord
    start <- sub_anno[i,]$UPA_start
    end <- sub_anno[i,]$UPA_end
    ftr <- sub_anno[i,]$ftr
    PA_id <- sub_anno[i,]$PA_id
    ftr_start <- sub_anno[i,]$ftr_start
    ftr_end <- sub_anno[i,]$ftr_end
    
    sub_db <- annoDB[annoDB$gene == gene_id,]
    ref_start <- sub_db$UPA_start
    ref_end <- sub_db$UPA_end
    
    
    #判断是否存在区间交集
    #obs=10表示一个范围误差
    #排序
    ranks <- rank(ref_start)
    ref_start <- ref_start[ranks] -obs
    ref_end <- ref_end[ranks] +obs
    #loc <- start < ref_end & end > ref_start
    loc <- ref_start < end & ref_end > start
    if(sum(loc) == 1){
      #当pa位点的起始和结束区间与参考坐标系的某些位点只存在一个重合区间时，即有唯一值
      sub_anno[i,]$PA_id2 <- sub_db[ranks,]$PA_id[loc]
    }else if(sum(loc) > 1){
      #当有多个重合区间时，通过距离计算筛选出唯一值，用参考坐标系中距离最近的pA位点进行注释
      ref_coord <- as.numeric(sub_db[ranks,]$coord[loc])
      loc2 <- which.min(abs(ref_coord-coord))
      sub_anno[i,]$PA_id2 <- sub_db[ranks,]$PA_id[loc][loc2]
    }else{
      dt <- sub_db[1,]
      dt$PA_id <- PA_id
      dt$start <- start
      dt$end <- end
      dt$coord <- coord
      dt$ftr_start <- ftr_start
      dt$ftr_end <- ftr_end
      dt$source <- source_name
      annoDB <- rbind(annoDB,dt)
    }
  }
  
  #返回原始的anno
  locs <- match(rownames(sub_anno),rownames(pa_anno))
  pa_anno$PA_id2[locs] <- sub_anno$PA_id2
  
  #DB整合进不包含的PA
  
  No_anno <- pa_anno[!(pa_anno$gene %in% annoDB$gene_id),]
  new_anno <- No_anno[,c("PA_id","chr","start","end","strand","coord","ftr","gene","ftr_start","ftr_end")]
  
  if(nrow(new_anno) >1){
    new_anno$source <- source_name
    colnames(new_anno)<- colnames(annoDB)
    annoDB <- rbind(annoDB,new_anno)
  }
  
  
  ##
  end_time <- Sys.time() 
  close(pb)
  run_time <- end_time - star_time
  
  return(list(pa_anno,annoDB))
  
}

##合并pA位点--------------------------------------------------------------------
##input data must contain "start”，“end","gene"
##example:
#PAC_id <- combind_PAC(utr3_anno)
#utr3_anno$PAC_id <- PAC_id

combind_PAC <- function(annos){
  PAC_id <- c()
  #计算数据框annos中gene列里不同基因的数量，这个数量用于后面设置进度条
  len <- length(table(factor(annos$gene)))
  pb <- txtProgressBar(style=3)  #ext progress bar in the R console，
  star_time <- Sys.time()  
  
  #开始一个循环，遍历annos数据中gene列的所有不同基因名称
  for(i in levels(factor(annos$gene))){
    #找到当前基因名称i在annos数据中的行索引
    dip <- which(levels(factor(annos$gene)) == i)
    setTxtProgressBar(pb, dip/len)   #根据当前处理的基因数目与总基因数目的比例，更新控制台中的进度条
    
    #提取当前基因i的注释信息
    tmp_anno <- annos[annos$gene == i,]
    ns <- tmp_anno[1,]$start #将当前基因的第一个pA注释的起始位置赋值给变量ns
    ne <- tmp_anno[1,]$end  #将当前基因的第一个pA注释的终止位置赋值给变量ne
    kk <- 1                 #初始化一个计数变量kk
    n_id <- paste0("PAC",i,"_",kk)  #生成一个基因的初始合并标识符
    PAC_id <- c(PAC_id,n_id)  #将初始合并标识符添加到PAC_id向量中
    
    #- 遍历当前基因的所有pA注释的位置区间，检查是否有重叠。
    #- 如果有重叠，则将相同的合并标识符n_id添加到PAC_id向量中，并根据需要更新ne值(即当同一个基因下一个pA结束位置大于上一个pA)。
    #- 如果没有重叠，则增加kk的值，并生成新的合并标识符n_id，更新ne和ns值。
    
    if(nrow(tmp_anno) > 1){
      for (j in 2:nrow(tmp_anno)) {
        if(tmp_anno[j,]$start <= ne){
          #overlap
          #set PAC_id
          n_id <- paste0("PAC",i,"_",kk)
          PAC_id <- c(PAC_id,n_id)
          #reset:ne
          if(tmp_anno[j,]$end > ne){
            ne <- tmp_anno[j,]$end 
          }
        }else{
          # no overlap
          #set PAC_id
          kk <- kk+1
          n_id <- paste0("PAC",i,"_",kk)
          PAC_id <- c(PAC_id,n_id)
          #reset:ne,ns
          ne <- tmp_anno[j,]$end 
          ns <- tmp_anno[j,]$start 
        }
      }
    }
  }
  end_time <- Sys.time() 
  close(pb)
  run_time <- end_time - star_time
  
  return(PAC_id)
}


###定义每个基因不重叠的3UTR区域--------------------------------------------------
##input data must contain "start","end","gene_id"
##example:
##utr3_bed=as.data.table(utr3_bed)
##setkey(utr3_bed,gene_id,start)
##new_utr=combind_3UTR(utr3_bed)

combind_3UTR <- function(annos){
  ##注释UTR区域，将区域有交叠的UTR区域注释为1utr_id
  print("Anno_unoverlap_utr")
  utr_id <- c()
  #计算数据框annos中gene列里不同基因的数量，这个数量用于后面设置进度条
  len <- length(unique(annos$gene_id))
  pb <- txtProgressBar(style=3)  #ext progress bar in the R console，
  star_time <- Sys.time()  
  
  #开始一个循环，遍历annos数据中gene列的所有不同基因名称
  for(i in 1:len){
    #找到当前基因名称i在annos数据中的行索引
    #dip <- which(levels(factor(annos$gene)) == i)
    gene=unique(annos$gene_id)[i]
    setTxtProgressBar(pb, i/len)   #根据当前处理的基因数目与总基因数目的比例，更新控制台中的进度条
    
    #提取当前基因的3UTR注释信息
    tmp_anno <- annos[annos$gene_id == gene,]
    ns <- tmp_anno[1,]$start #将当前基因的第一个3UTR注释的起始位置赋值给变量ns
    ne <- tmp_anno[1,]$end  #将当前基因的第一个3UTR注释的终止位置赋值给变量ne
    kk <- 1                 #初始化一个计数变量kk
    n_id <- paste0(gene,"_",kk)  #生成一个基因的初始合并标识符
    utr_id <- c(utr_id,n_id)  #将初始合并标识符添加到utr_id向量中
    
    #- 遍历当前基因的所有UTR注释的位置区间，检查是否有重叠。
    #- 如果有重叠，则将相同的合并标识符n_id添加到utr_id向量中，并根据需要更新ne值(即当同一个基因下一个utr结束位置大于上一个转录本的utr)。
    #- 如果没有重叠，则增加kk的值，并生成新的合并标识符n_id，更新ne和ns值。
    
    if(nrow(tmp_anno) > 1){
      for (j in 2:nrow(tmp_anno)) {
        if(tmp_anno[j,]$start <= ne){
          #overlap
          #set PAC_id
          n_id <- paste0(gene,"_",kk)
          utr_id <- c(utr_id,n_id)
          #reset:ne
          if(tmp_anno[j,]$end > ne){
            ne <- tmp_anno[j,]$end 
          }
        }else{
          # no overlap
          #set PAC_id
          kk <- kk+1
          n_id <- paste0(gene,"_",kk)
          utr_id <- c(utr_id,n_id)
          #reset:ne,ns
          ne <- tmp_anno[j,]$end 
          ns <- tmp_anno[j,]$start 
        }
      }
    }
  }
  end_time <- Sys.time() 
  close(pb)
  run_time <- end_time - star_time
  print(run_time)
  
  ##根据utr_id,合并区域有交叠的utr区域
  print("Merge_overlap_utr")
  annos$utr_id=utr_id
  
  unique_id=unique(utr_id)
  new_utr=data.frame()
  pb2 <- txtProgressBar(style=3)  #ext progress bar in the R console，
  star_time2 <- Sys.time() 
  
  for(k in 1:length(unique_id)){
    setTxtProgressBar(pb2, k/length(unique_id))
    
    id=unique_id[k]
    temp=subset(annos,annos$utr_id==id)
    temp_utr=data.frame(chromosome=unique(temp$chromosome),
                        start=min(temp$start),
                        end=max(temp$end),
                        gene_id=unique(temp$gene_id),
                        utr_id=id,
                        strand=unique(temp$strand))
    new_utr=rbind(new_utr,temp_utr)
    
  }
  end_time2 <- Sys.time() 
  close(pb2)
  run_time2 <- end_time2 - star_time2
  print(run_time2)
  
  
  return(new_utr)
}


#计算pA位点与参考pA位点之间最近的距离-----------------------------------
getDistance.region  <- function(query.pa=NULL,known.pa=NULL,srrID=NULL,group=NULL,species=NULL,protocol=NULL){
  distance.pa <- GenomicRanges::distanceToNearest(query.pa,known.pa)
  #length(distance.pa );length(query.pa);length(known.pa);length(unique(distance.pa@from ))
  
  temp.data <- data.frame(distance=distance.pa@elementMetadata$distance,
                          group=group,srrID,
                          species=species,
                          protocol=protocol)
  #index1：计算 query.pa 中每个元素的结束位置是否在 known.pa 对应元素的结束位置之前
  index1 <- which(end(query.pa)[distance.pa@from]-end(known.pa)[distance.pa@to]<0)
  #index2：找出 query.pa 中链方向为 "+" 的元素的索引
  index2 <- which(as.character(strand(query.pa)[distance.pa@from]) =="+")
  #index12：找出同时满足上述两个条件的元素的索引
  index12 <- intersect(index1,index2)
  
  #index3：计算 query.pa 中每个元素的结束位置是否在 known.pa 对应元素的结束位置之后
  index3 <- which(end(query.pa)[distance.pa@from]-end(known.pa)[distance.pa@to]>0)
  #index4：找出 query.pa 中链方向为 "-" 的元素的索引
  index4 <- which(as.character(strand(query.pa)[distance.pa@from]) =="-")
  #index34：找出同时满足上述两个条件的元素的索引
  index34<- intersect(index3,index4)
  #temp.data$distance_raw <- temp.data$distance
  #调整距离：对于在+/-链在位于known.pa 对应元素的结束位置之前/后的距离调整为负值，表示位于known.pa的上游位置
  temp.data$distance[c(index12,index34)] <- temp.data$distance[c(index12,index34)]*(-1)
  
  #有时查询区间在所有已知区间中没有找到足够近的区间，因此没有返回距离值，故先删除没有查询到与refpA的最近距离值的pA位置信息
  query.pa=query.pa[distance.pa@from]
  
  temp.data$count <- query.pa$count
  query.pa$distance <- temp.data$distance
  return(list(distance.data=temp.data,query.data=query.pa))
  
}


#绘制柱状图表示重叠的pA在每种工具中所占的比例-----------------------------------
plotNovel <- function(novel.data=NULL,input.level=NULL){
  novel.data$pro <- round(novel.data$Freq/novel.data$total,4)
  
  
  novel.data$type<- factor(novel.data$type,
                           levels = c("Novel","Overlapping"))
  
  novel.data <- novel.data[order(novel.data$type,decreasing = T),]
  #对于每个工具，计算Freq列的累积和，并将结果赋值给新列label_ypos
  novel.data <- plyr::ddply(novel.data, "group",
                            transform, label_ypos=cumsum(Freq))
  #novel.data$group <- factor(novel.data$group,
  #                           levels=input.level)
  
  
  
  novel.data$pro <- paste0(round(novel.data$pro*100,2),"%")
  novel.data$pro[novel.data$type=="Novel"] <- ""
  
  #输出每种工具识别的pA中重叠的pA的数量和比例，并使其比例显示在柱状图上
  novel.data_sub<- subset(novel.data,novel.data$type=="Overlapping")
  #novel.data_sub$group <- factor(novel.data_sub$group,
  #                               levels=input.level)
  
  
  p<- ggplot(data=novel.data,
             aes(x=group, y=Freq,fill=type)) +
    geom_bar(stat="identity")+p_theme+
    # fill=rep(c("#338DBF","#0071AF"),time=length(unique(novel.data$source))))+p_theme+
    labs(x=NULL,y=input.level)+
    guides(fill=guide_legend(title=NULL))+#scale_fill_brewer(palette="Set3")+
    geom_text(data=novel.data_sub,aes(y=Freq, label=pro), vjust=2.6, 
              color="white", size=3)+
    theme(legend.position = "top",axis.text.x = element_text(angle = 45,hjust=1))+
    scale_fill_manual(values=rev(wes_palette(n=2, name="Royal1")))
  return(list(p=p,novel.data=novel.data))
}

##绘制热图气泡图展示不同样本数据中具有refpA支持的比例-------------------------------
#input.label:梯度热图标签，表示每种工具识别的真实pA(与refpa之间的距离小于50bp)的比例(%)
#input.label2：圆点大小标签，表示每种工具识别的总的3UTR pA的数目(x10^3)
plot.Heatmap <- function(pac.data.need = NULL,input.label = "real 3UTR PAC%",
                         input.label2 =bquote("# 3UTR PAC (x 10"^3~")"),
                         sample.info=NULL){
  
  p<- ggplot(pac.data.need, aes(y = factor(srr),
                                x = factor(group))) +        ## global aes
    #geom_tile(aes(fill = rectheat),color="black",fill="#F5F5F5") + 
    #      ## to get the rect filled
    geom_point(aes(fill = pro.value, 
                   size =total),  shape = 21,colour= "#C7B19C",stroke = 1.2)  +    ## geom_point for circle illusion;shape=21:实心圆形；shape=23,实心菱形
    scale_fill_gradientn(
      colours = heatmap.col4,
      name=input.label)+
    #scale_fill_gradientn(colours=c("#023FA5","#7D87B9","#BEC1D4",
    #                                  rev(colorspace::heat_hcl(20))),name=input.label)+       ## color of the corresponding aes
    #scale_size(range = c(1, 20))+             ## to tune the size of circles
    theme_classic()+
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
      check_overlap=T,color="black",
      hjust=-0.7#vjust=-1
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
  
  
  p1.special <- ggplot(temp, aes(y=srr, x=anno2, fill=species)) + geom_tile() + 
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

#######################################################################
plot.Heatmap2 <- function(pac.data.need = NULL,input.label = "real 3UTR PAC%",
                         input.label2 =bquote("# 3UTR PAC (x 10"^3~")"),
                         sample.info=NULL){
  
  p<- ggplot(pac.data.need, aes(y = factor(srr),
                                x = factor(group))) +        ## global aes
    geom_tile(aes(fill = rectheat),color="black",fill="#F5F5F5") + 
    #      ## to get the rect filled
    geom_point(aes(fill = pro.value, 
                   size =total),  shape = 21,colour= "#C7B19C",stroke = 1.2)  +    ## geom_point for circle illusion;shape=21:实心圆形；shape=23,实心菱形
    scale_fill_gradientn(
      colours = heatmap.col4,
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
  
  
  p1.special <- ggplot(temp, aes(y=srr, x=anno2, fill=species)) + geom_tile() + 
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

###########################################################################


#统计基因组不同区域的pA数量并计算poly(A)信号的分布----------------------------------
##结果返回列表，包含以下三个元素：1.data.PACds（注释的pa位点，包含每个pA与临近refpA之间的距离）
##2.feature.stastic：pA在基因组不同区域的分布（pro_feature列：基因组不同区域的pA所占的比例/
##pro_type列：在基因组不同区域，novel pa与overlap pa在所有novel pa/overlap pa种所占的比例）
##3.poly.signal：polyA信号的·分布

##get.annotation:注释pA，统计基因组不同区域的pA数量并计算poly(A)信号的分布
##参数species：Human/Mouse/Arabidopsis,当species=Human/Mouse/Arabidopsis且ext3UTRlen=NULL,则默认ext3UTRlen=2000/2000/1000；
#若ext3UTRlen不为NULL,则ext3UTRlen未指定值
get.annotation<- function(input.data=NULL,
                          annotation=NULL,
                          species=NULL,
                          ext3UTRlen=NULL,
                          source=NULL,
                          cutoff.distance=50,
                          bsgenome=NULL){
  input.data <- data.frame(input.data)
  #head(input.data )
  colnames(input.data )[1:3] <- c("chr","UPA_start","UPA_end")
  input.data$type_distance <- "Novel"
  input.data$type_distance[which(abs(input.data$distance)<=cutoff.distance)] <- "Overlapping"
  input.data=select(input.data,-gene)
  
  data.PACds <- readPACds(input.data, colDataFile=NULL)
  data.PACds <- annotatePAC(data.PACds, aGFF = annotation)
  if(species=="Human"& is.null(ext3UTRlen)){
    data.PACds <- ext3UTRPACds(data.PACds,ext3UTRlen = 2000)
  }
  if(species=="Mouse"& is.null(ext3UTRlen)){
    data.PACds <- ext3UTRPACds(data.PACds,ext3UTRlen = 2000)
  }
  if(species=="Arabidopsis"& is.null(ext3UTRlen)){
    data.PACds <- ext3UTRPACds(data.PACds,ext3UTRlen = 1000)
    #修改染色体名称，与bsgenome对象中的染色体名称匹配
    data.PACds@anno$chr=paste0("Chr",data.PACds@anno$chr)
  }
  if(!is.null(ext3UTRlen)){
    data.PACds <- ext3UTRPACds(data.PACds,ext3UTRlen = ext3UTRlen)
  }
  
  data.PACds<- annotateByPAS(data.PACds, bsgenome, grams='AATAAA', from=-50, to=1, label=NULL)
  data.PACds <- annotateByPAS(data.PACds, bsgenome, grams='V1', from=-50, to=1, label=NULL)
  data.PACds@anno$pA.signal <- "Others"
  data.PACds@anno$pA.signal[which(!is.na(data.PACds@anno$V1_dist))] <- "1Variants"
  data.PACds@anno$pA.signal[which(!is.na(data.PACds@anno$AATAAA_dist))] <- "AATAAA"
  
  #基因组不同区域的novel/overlap pa附近的polyA信号分布
  pa.label <- as.data.frame( table(paste0(data.PACds@anno$pA.signal,
                                          "_",data.PACds@anno$ftr,"_",data.PACds@anno$type_distance) ))
  colnames(pa.label) <- c("Features1","No")
  
  pa.label$pa.signal <- as.character(limma::strsplit2(pa.label$Features1,"\\_")[,1])
  pa.label$Features <- as.character(limma::strsplit2(pa.label$Features1,"\\_")[,2])
  pa.label$type <- as.character(limma::strsplit2(pa.label$Features1,"\\_")[,3])
  pa.label$source <- source
  
  #table(data.PACds@anno$pA.signal)
  #1Variants    AATAAA    Others 
  #37644     56329      8970 
  data.stastic <- as.data.frame( table(paste0(data.PACds@anno$ftr,"_",data.PACds@anno$type_distance) ))
  colnames(data.stastic) <- c("Features1","No")
  data.stastic$Features <- as.character(limma::strsplit2(data.stastic$Features1,"\\_")[,1])
  data.stastic$type <- as.character(limma::strsplit2(data.stastic$Features1,"\\_")[,2])
  
  label.need <- c("3UTR","exon",
                  "intron","5UTR","CDS","intergenic")
  label.need1 <- c("3UTR","exon",
                   "intron","5UTR","exon","intergenic")
  data.stastic$Features <- factor(data.stastic$Features,
                                  levels=c(label.need,setdiff(data.stastic$Features,label.need)),
                                  labels=c(label.need1,setdiff(data.stastic$Features,label.need)))
  
  #sum_novel total #novel and overlapping pA
  data.stastic <- data.stastic %>%
    dplyr::group_by(type) %>% 
    dplyr::mutate(sum_novl=sum(No))
  #sum_features total #features 
  data.stastic <- data.stastic %>%
    dplyr::group_by(Features) %>% 
    dplyr::mutate(sum_feature=sum(No))
  #pro_feature:基因组不同区域的pA所占的比例
  data.stastic$pro_feature <- round(data.stastic$sum_feature/sum(data.stastic$No),4)
  #pro_type：在基因组不同区域，novel pa与overlap pa在所有novel pa/overlap pa种所占的比例
  data.stastic$pro_type <- round(data.stastic$No/data.stastic$sum_novl,4)
  data.stastic$source <- source
  return(list(data.PACds=data.PACds,feature.stastic=data.stastic,poly.signal=pa.label))
}

############################################################
#统计基因组不同区域的pA数量
data.stastic.temp <- function(input.data=NULL){
  data.stastic <- as.data.frame( table(paste0(input.data$ftr,"_",input.data$type_distance) ))
  colnames(data.stastic) <- c("Features1","No")
  data.stastic$Features <- as.character(limma::strsplit2(data.stastic$Features1,"\\_")[,1])
  data.stastic$type <- as.character(limma::strsplit2(data.stastic$Features1,"\\_")[,2])
  label.need <- c("3UTR","exon",
                  "intron","5UTR","CDS","intergenic")
  label.need1 <- c("3UTR","exon",
                   "intron","5UTR","exon","intergenic")
  data.stastic$Features <- factor(data.stastic$Features,
                                  levels=c(label.need,setdiff(data.stastic$Features,label.need)),
                                  labels=c(label.need1,setdiff(data.stastic$Features,label.need)))
  
  data.stastic <- data.stastic %>%
    dplyr::group_by(type) %>% 
    dplyr::mutate(sum_novl=sum(No))
  #sum_features total #features 
  data.stastic <- data.stastic %>%
    dplyr::group_by(Features) %>% 
    dplyr::mutate(sum_feature=sum(No))
  
  data.stastic$pro_feature <- round(data.stastic$sum_feature/sum(data.stastic$No),4)
  data.stastic$pro_type <- round(data.stastic$No/data.stastic$sum_novl,4)
  return(data.stastic)
}


##绘制单核苷酸分布图------------------------------------------
plotATCGforFAfile.cp <- function (faFiles, ofreq=FALSE, opdf=TRUE, refPos=NULL, start=NA, end=NA, filepre='', mergePlots=FALSE) {
  
  #library(Biostrings, verbose = FALSE)
  #library(ggplot2)
  
  if (filepre=='' & (opdf & mergePlots)) stop("mergePlots & opdf=TRUE, should provide output file name: filepre\n")
  
  dmerge=c()
  
  if (opdf) {
    if (mergePlots) {
      pdf(file=paste0(filepre,'.pdf'))
    }
  }
  
  
  for (fafile in faFiles) {
    seq=readDNAStringSet(fafile, format="fasta")
    nseq=length(seq)
    if (!is.na(start) | !is.na(end)) {
      seq=subseq(seq, start=start, end=end)
      refPos=NULL
    }
    din<-consensusMatrix(seq, as.prob=TRUE, baseOnly=TRUE)
    din<- as.data.frame(t(din))
    din$other <- NULL
    din$pos <- c(1:nrow(din))
    if (!is.null(refPos)) {
      din$pos=din$pos-refPos
    }
    opre=gsub('\\.fa$|\\.fasta$','',fafile, ignore.case = TRUE)
    omain=paste0(opre, ' #',nseq)
    if (filepre!='') {
      opre=paste0(opre,'.',filepre)
    }
    if (ofreq) {
      din=din[, c('pos','A','C','G','T')]
      ofreqfile=paste0(opre,'.freq')
      write.table(din, file=ofreqfile, col.names = TRUE, row.names = FALSE, sep="\t", quote=F)
      cat('>>>',ofreqfile,'\n')
    }
    
    if (mergePlots) {
      if (length(dmerge)==0) {
        dmerge=cbind(source=omain, din)
      } else {
        dmerge=rbind(dmerge, cbind(source=omain,din))
      }
    }
    din <- reshape2::melt(din,id.vars='pos',variable.name = 'base', value.name = 'freq')
    ATCGpicture<-ggplot(din,aes(x=din$pos,y=din$freq,colour=din$base)) + geom_line() +
      xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw() +
      ggtitle(omain)
    
    if (opdf & !mergePlots) {
      pdf(file=paste0(opre,'.pdf'))
      return(ATCGpicture)
      dev.off()
      cat('>>>',paste0(opre,'.pdf'),'\n')
    }
    ## only plot the last one
    if (!opdf & !mergePlots & faFiles[length(faFiles)]==fafile) print(ATCGpicture)
  }
  
  if (mergePlots) {
    din <- reshape2::melt(dmerge, id.vars=c('source','pos'), variable.name = 'base', value.name = 'freq')
    din$source=as.character(din$source)
    #if file name is too long, remove common chars
    #if no common chars, change the first 1-30 chars to ...
    l=min(nchar(din$source))
    # if (l>=30) {
    #   din$source =.autoDetectCommonString(unique(din$source), sbj=din$source, beAll=TRUE)
    #   l=min(nchar(din$source))
    #   if (l>=30) {
    #     din$source=paste0('...', substr(din$source, 31, nchar(din$source)))
    #   }
    # }
    
    ATCGpicture <- ggplot(din, aes(x=pos, y=freq, group=base)) +
      geom_line(aes(colour=base)) +
      xlab("Position") + ylab('Base component') + theme(legend.title=element_blank()) + theme_bw()
    facet_grid(. ~ source)
    ATCGpicture = ATCGpicture + facet_wrap(~ source, ncol=3)
    
    return(list(ATCGpicture=ATCGpicture,din=din))
    if (opdf) {
      dev.off()
      cat('>>>',paste0(filepre,'.pdf'),'\n')
    }
  }
  
}

########plot.UMAP------------------------------------------------------------------
plotUMAP <- function(data,group=NULL,allcolour=NULL,title=NULL,xend_add=2,yend_add=4,annotate=TRUE,label_size=2,point.size=1,label=TRUE){
  library(ggrepel)
  library(conflicted)
  library(ggsci)
  if(class(data) == "Seurat"){
    umap <- data@reductions$umap@cell.embeddings %>%  #坐标信息
      as.data.frame() %>% 
      cbind(cell_type = data@meta.data[,group]) # 注释后的label信息 ，改为cell_type
  }else if(class(data) == "data.frame"){
    umap <- data
  }
  
  if(is.null(allcolour)){
    #allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
    p <- ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  
      geom_point(size = 1 , alpha =1 )  + 
      scale_color_d3(palette = "category20")+
      labs(title=title)
    
  }else{
    p=ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
      geom_point(size = point.size, alpha = 1) +
      scale_color_manual(values = allcolour) +  # 使用自定义颜色
      labs(title = title)  
    }
  
  
  
  # 去掉网格线，坐标轴和背景色即可
  p2 <- p  +
    theme(panel.grid.major = element_blank(), #主网格线
          panel.grid.minor = element_blank(), #次网格线
          panel.border = element_blank(), #边框
          axis.title = element_blank(),  #轴标题
          axis.text = element_blank(), # 文本
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white'), #背景色
          plot.background=element_rect(fill="white"))
  
  
  # legeng部分去掉legend.title后，调整标签大小，标签点的大小以及 标签之间的距离
  p3 <- p2 +         
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(), #去掉legend.title 
      legend.key=element_rect(fill='white'), #
      legend.text = element_text(size=7), #设置legend标签的大小
      legend.key.size=unit(0.5,'cm') ) +  # 设置legend标签之间的大小
    guides(color = guide_legend(override.aes = list(size=4))) #设置legend中 点的大小 
  
  
  # 坐标轴放到左下角可以通过ggplot2添加箭头和文本实现。
  
  conflict_prefer("geom_segment", "ggplot2")
  if(annotate){
    p4 <- p3 + geom_segment(aes(x = min(UMAP_1),y = min(UMAP_2),
                                xend = min(UMAP_1)+xend_add, yend = min(UMAP_2)),
                            colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
      geom_segment(aes(x = min(UMAP_1)  , y = min(UMAP_2),
                       xend = min(UMAP_1) , yend = min(UMAP_2) + yend_add),
                   colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
      annotate("text", x = min(umap$UMAP_1) +2.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
               color="black",size = 4, fontface="bold" ) + 
      annotate("text", x = min(umap$UMAP_1) -1.5, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
               color="black",size = 4, fontface="bold" ,angle=90) 
  }else{
    p4<-p3
  }
  # 添加 label
  
  # 调整umap图 - repel - labels
  
  cell_type_med <- umap %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  p5 <- p4 +geom_label_repel(aes(label=cell_type),
                             data = cell_type_med,
                             size=label_size, fontface="bold",
                             box.padding=unit(0.5, "lines"), point.padding=unit(1.6, "lines"),
                             segment.size = 0.5,
                             arrow = arrow(length=unit(0.01, "npc")),force = 1, max.iter = 3e3) +
    theme(legend.position = "none")
  if(label==TRUE){
    return(p5)
  }else{
    return(p4)
  }
  
}


##calculate_wf-----------------------------------------
# 加权频次（WF）计算函数
# 输入：deg_lists - 包含各工具DE APA gene以及显著性排名和百分比排名的列表
# 输出：数据框，包含基因名称和对应的WF值
calculate_wf<- function(deg_lists) {
  # 获取所有唯一基因
  all_genes <- unique(unlist(lapply(deg_lists, function(df) df[, 1])))
  
  # 创建结果数据框
  wf_df <- data.frame(
    gene = all_genes,
    wf = numeric(length(all_genes)),
    stringsAsFactors = FALSE
  )
  
  # 计算每个基因的WF值
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]
    wf_score <- 0
    
    # 遍历每个工具
    for (tool_name in names(deg_lists)) {
      # 获取该工具的DE APA gene列表和完整排名数据
      tool_deg <- deg_lists[[tool_name]]
      
      # 检查基因是否在该工具的DE APA gene列表中
      if (gene %in% tool_deg$gene_id) {
        
        # 获取基因在该工具中的百分位数排名PR排名
        pr<- tool_deg[tool_deg$gene_id == gene,]$PR
        
        
        # 计算贡献值（1 - PR/100转换为比例）
        contribution <- 1 - (pr / 100)
        
        wf_score <- wf_score + contribution
      }
    }
    
    #记录每个基因最后的wf得分
    wf_df$wf[i] <- wf_score
  }
  
  # 按WF值降序排序
  wf_df <- wf_df[order(-wf_df$wf), ]
  
  return(wf_df)
}





###****Difference Analysis***********------------------------------------------------------

  
  PreInput <- function(counts, samples, tx2gene) {
    list(
      counts = counts,
      samples = samples,
      tx2gene = tx2gene
    )
  }

# 预处理输入数据
#input_data <- PreInput(counts=cts, samples, tx2gene)


####run_DEXSeq-----------------------------------------------------------------------

run_DEXSeq<- function(input_data,conditions){
  
  input_data <- list(
    counts = input_data$counts,
    samples = input_data$samples,
    tx2gene = input_data$tx2gene
  )
  
  rownames(input_data$counts) <- input_data$counts$tx_id
  counts <- input_data$counts[, -(1:2)] 
  
  samples <- data.frame(
    sample_id = input_data$samples$sample_id,
    group = input_data$samples$condition,
    stringsAsFactors = FALSE
  )
  
  #colnames(input_data$counts)[colnames(input_data$counts) == "tx_id"] <- "feature_id"
  #dmdata <- DRIMSeq::dmDSdata(counts = input_data$counts, samples = input_data$samples)
  
  #构建 DEXSeqDataSet 对象
  #colData:包含样本分组信息
  #counts:表达矩阵
  #rowRanges:转录本的注释信息，需包含转录本（featureID/exonID)对应的基因信息（groupID）
  dxd <- DEXSeq::DEXSeqDataSet(
    countData = round(counts,0),  #要求原始的表达矩阵计数，不适合标准化后的计数,这里先四舍五入取整
    sampleData = samples,
    design = ~sample + exon + group:exon, #用于建模差异外显子使用（Differential Exon Usage, DEU）分析的统计模型
    #~sample:控制样本间的基线差异（如测序深度、文库大小等）
    #/exon:控制不同外显子的基线表达差异/group:exon 检测外显子使用比例在不同组间的差异（核心目标）
    featureID = input_data$counts$tx_id,
    groupID = input_data$counts$gene_id
  )
  
  ##这里使用原本的TPM值重新储存在counts对象中
  #DEXSeq 的差异分析基于 负二项分布，要求输入为原始计数或缩放计数，而TPM是连续值，破坏了模型的离散性假设。
  #使用TPM会导致假设检验（如 testForDEU）失效，可能产生假阳性/阴性
  
  #head(counts(dxd))  #注：这里的列数会变为原本的两倍，前一半的列为原本转录本的表达计数，后一半为同一基因上其他转录本的表达计数总和
  #row.names(counts)=paste0(input_data$counts$gene_id,":",input_data$counts$tx_id)
  #counts=counts[row.names(counts(dxd)),]
  #assay(dxd, withDimnames=FALSE) <- as.matrix(counts)
  
  #标准化：校正样本间的 测序深度差异(检测外显子使用差异（exon usage）)
  dxd <- DEXSeq::estimateSizeFactors(dxd)
  #head(colData(dxd)),sizeFactor这一列为标准化后的值
  
  ##估计离散度
  
  dxd <- tryCatch({
    # 第一次尝试：默认参数（fitType="parametric"）
    DEXSeq::estimateDispersions(dxd)
  }, error = function(e) {
    message("Default parametric fit failed (", e$message, "), trying local fit...")
    
    # 第二次尝试：局部回归（fitType="local"）
    tryCatch({
      DEXSeq::estimateDispersions(dxd, fitType = "local", quiet = TRUE)
    }, error = function(e) {
      message("Local fit failed (", e$message, "), switching to mean fit...")
      
      # 第三次尝试：均值拟合（fitType="mean"）
      DEXSeq::estimateDispersions(dxd, fitType = "mean", quiet = TRUE)
    })
  })
  
  
  
  head(rowData(dxd))
  #plotDispEsts(dxd)  # 绘制离散度与表达量均值的关系
  #黑点：每个外显子的原始离散度估计。/红线：拟合的离散度趋势（均值-离散度关系）。
  #最终离散度：根据趋势调整后的离散度值（用于后续检验）。
  
  #(2) 典型模式
  #高表达外显子：离散度较低（技术噪音占比小，变异主要由生物学因素驱动）。
  #低表达外显子：离散度较高（技术噪音主导，数据更分散）。
  
  #(3) 异常情况
  #离散度过高：可能表明样本存在批次效应或技术问题。/拟合线不平稳：需检查数据质量或模型假设。
  
  #移除 condition:exon 交互项，即 忽略实验条件对外显子使用的影响。
  #若完整模型与简化模型拟合差异显著，则说明条件对外显子使用有影响（存在差异外显子使用）
  
  dxd <- DEXSeq::testForDEU(dxd, reducedModel = ~ sample + exon)
  #估计外显子在不同条件间使用率变化
  #colData(dxd)
  #指定对比关系,DEXSeq会将因子水平在前的指定为对照组，故groups[2]在前，为对照组
  colData(dxd)$group <- factor(colData(dxd)$group, 
                               levels = c(conditions[2],conditions[1]))
  dxd = DEXSeq::estimateExonFoldChanges( dxd, fitExpToVar="group")
  
  # 提取结果
  dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)
  ResDEXSeq<-as.data.frame(dxr)
  colnames(ResDEXSeq)[colnames(ResDEXSeq) == "featureID"] <- "tx_id"
  colnames(ResDEXSeq)[colnames(ResDEXSeq) == "groupID"] <- "gene_id"
  
  ##  [1] "group/gene identifier"                                       
  ##  [2] "feature/exon identifier"                                     
  ##  [3] "mean of the counts across samples in each feature/exon"      
  ##  [4] "exon dispersion estimate"                                    
  ##  [5] "LRT statistic: full vs reduced"                              
  ##  [6] "LRT p-value: full vs reduced"                                
  ##  [7] "BH adjusted p-values"                                        
  ##  [8] "exon usage coefficient"                                      
  ##  [9] "exon usage coefficient"                                      
  ## [10] "relative exon usage fold change"                             
  ## [11] "GRanges object of the coordinates of the exon/feature"       
  ## [12] "matrix of integer counts, of each column containing a sample"
  
  #head(ResDEXSeq)
  #保留前10列的结果
  ResDEXSeq=ResDEXSeq[,c(1:10)]
  
  
  
  #返回结果
  return(ResDEXSeq)
}



###构建伪批量函数--------------------------------------------------------------------
create_pseudo_bulk<- function(counts, 
                              samples,
                              population.1 = NULL, 
                              population.2 = NULL, 
                              num.splits = 6, 
                              seed.use = 1,
                              replicates.1 = NULL,
                              replicates.2 = NULL){
  set.seed(seed.use)
  apa.seurat.object <- Seurat::CreateSeuratObject(counts = counts,meta.data = samples)
  peaks.use = rownames(apa.seurat.object)
  ## 设置细胞标识
  
  Seurat::Idents(apa.seurat.object) <- samples$condition
  if (is.null(replicates.1)) {
    
    if (length(population.1) == 1) {
      cells.1 <- names(Seurat::Idents(apa.seurat.object))[which(Seurat::Idents(apa.seurat.object) == population.1)]
    } else{
      cells.1 <- population.1
    }
    
    cells.1 = sample(cells.1)   # 随机打乱细胞顺序
    cell.sets1 <- split(cells.1, sort(1:length(cells.1)%%num.splits)) #均分为6个样本组
  } else{
    ## user has provided cells for replicates - use these instead
    cell.sets1 <- replicates.1
  }
  
  ## create a profile set for first cluster
  profile.set1 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- Seurat::GetAssayData(apa.seurat.object, slot = "counts", assay = "RNA")[peaks.use, this.set] #提取目标转录本（peaks.use）和目标细胞（this.set）
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x))) # 多细胞时求和
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix #单个细胞时直接赋值
    }
  }
  rownames(profile.set1) <- peaks.use
  colnames(profile.set1) <- paste0("Population1_", 1:length(cell.sets1))
  
  ## create a profile set for second cluster
  if (is.null(replicates.2)) {
    if (is.null(population.2)) {
      cells.2 <- setdiff(colnames(apa.seurat.object), cells.1)
    } else {
      if (length(population.2) == 1) {
        #若直接给"clusterA"（单个簇名）	从 Idents 提取属于 clusterA 的细胞。
        cells.2 <- names(Seurat::Idents(apa.seurat.object))[which(Seurat::Idents(apa.seurat.object) == population.2)]
      } else {
        #若给细胞向量：c("cell1", "cell2")（向量），直接使用这些细胞
        cells.2 <- population.2
      }
    }
    
    cells.2 = sample(cells.2)
    cell.sets2 <- split(cells.2, sort(1:length(cells.2)%%num.splits))
  } else{
    ## user has provided cells for replicates - use these instead
    cell.sets2 <- replicates.2
  }
  
  
  profile.set2 = matrix(, nrow = length(peaks.use), ncol = length(cell.sets2))
  for (i in 1:length(cell.sets2)) {
    this.set <- cell.sets2[[i]]
    sub.matrix <- Seurat::GetAssayData(apa.seurat.object, slot = "counts", assay = "RNA")[peaks.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) sum(x)))
      profile.set2[, i] <- this.profile
    } else {
      profile.set2[, i] <- sub.matrix
    }
  }
  rownames(profile.set2) <- peaks.use
  colnames(profile.set2) <- paste0("Population2_", 1:length(cell.sets2))
  
  ## merge the count matrices together
  peak.matrix <- cbind(profile.set1, profile.set2)
  
  return(peak.matrix)
}








