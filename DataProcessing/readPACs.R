###****************read and annotation PACs*****************#####
##need R package
library(Sierra)
library(movAPA)
library(GenomicFeatures)
library(Seurat)
library(dplyr)
library(scAPA)
library(SCAPE)
library(ggplot2)
library(eoffice)
library(aqp)
library(Rsamtools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Athaliana.TAIR.TAIR9)
source("C:/R data/plot.function/plotATCGforFAfile.2.R")
source("G:/scAPA/code/function.R")
#source("G:/scAPA/scAPA_tools/scAPAtrap-master/R/scAPAtrap_funlib.R")

##note:When running different tools, be sure to modify the tool name and input file.------------------------

toolfile.list=c("scapatrap","sierra","scapture","infernape","scape","scapa","maaper","polyapipe","scraps","scutrquant")
tool.list=c("scAPAtrap","Sierra","SCAPTURE","Infernape","SCAPE","scAPA","MAAPER","polyApipe","scraps","scUTRquant")


##sample data---------------------------------------------------
GSElists=c("mouse_sperm","mouse_tcell","mouse_tip","mouse_esc","mouse_gfp","mouse_bone","mouse_hyp","mouse_int","human_pbmc","human_covid_pbmc",
           "arab")

datalists=c("merge","merge","merge","merge","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399","GSM2333586","GSM2333587",
            "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
            "GSM3490689","GSM3490690","GSM3490691")



##Data from different species, parameter settings-------
species="arab" #can be mouse/human/arab


if(species=="mouse"){
  reference.file <- "G:/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf"#mouse
  #gff<-parseGenomeAnnotation(reference.file)
  gff=makeTxDbFromGFF(reference.file , format="gtf")
  gff<-parseGenomeAnnotation(gff)
  
  
  chrs=paste0("chr",c(as.character(1:19),"X","Y"))#mouse
  bsgenome=BSgenome.Mmusculus.UCSC.mm10
  genome="mm10"
  ext.len=2000
  v=getVarGrams("mm")
}

if(species=="human"){
  reference.file <- "G:/ref_file/GRch38fa+gtf/v34/gencode.v34.annotation.gtf"
  gff=makeTxDbFromGFF(reference.file , format="gtf")
  gff<-parseGenomeAnnotation(gff)
  
  chrs=paste0("chr",c(as.character(1:22),"X","Y"))#human
  bsgenome=BSgenome.Hsapiens.UCSC.hg38
  genome="hg38"
  ext.len=2000
  v=c("AATAAA" ,"ATTAAA" ,"TATAAA", "AGTAAA" ,"AATACA" ,"CATAAA", "AATATA", "GATAAA" ,"AATGAA", "AATAAT" ,"AAGAAA" ,"ACTAAA" ,
      "AATAGA", "ATTACA", "AACAAA","ATTATA", "AACAAG", "AATAAG")
  
}

if(species=="arab"){
  reference.file <- "G:/ref_file/TAIR/Arabidopsis_thaliana.TAIR10.49.gtf"
  gff=makeTxDbFromGFF(reference.file , format="gtf")
  gff<-parseGenomeAnnotation(gff)
  
  chrs=c(as.character(1:5))#arab
  fafile<-"G:/ref_file/TAIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
  bsgenome=FaFile(fafile)
  indexFa(bsgenome$path)
  genome="tair10"
  ext.len=1000
  v=getVarGrams("V1")
  
}



GSE="arab"
GSEfile=paste0("G:/scAPA/PAC/",GSE)
datas=c("GSM3490689","GSM3490690","GSM3490691")


###process scAPAtrap PACs file-----------------------------------------


toolfile=toolfile.list[1]
tool=tool.list[1]

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)

file.path=paste0(GSEfile,"/",toolfile)

GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
#GSMlists="merge" #whenGSE=GSElists[1:4]
#GSMlists=datas

  for(i in 1: length(GSMlists)){
    message("start") 
    GSM=GSMlists[i]
    #path=paste0(file.path,"/",GSM,"/","APA.tails")#when tools is scapatrap
    path=paste0(file.path,"/",GSM,"/","APA.notails")
    setwd(path)
    
    
    if (!dir.exists("./PACds")){
      dir.create("./PACds")
    } else {
      print("Output dir already exists!")
    }
    
    
    ##filtered pa-----
    pafile=paste0(file.path,"/",GSM,"/","APA.tails/scAPAtrapData.rda")
    
    sub_scPACds=unify.pa(pafile = pafile,tool="scAPAtrap")
    
    ##raw pA-----------
    
    outputfile=paste0(file.path,"/",GSM,"/","APA.notails/scAPAtrapData.rda")
    
    
    sink("PACds/readPACs.log.txt")
    #将pa读入PACdataset对象
    scPACds=unify.pa(pafile = outputfile,tool="scAPAtrap")
    
   
    if(species=="arab"){
     scPACds@anno$chr=gsub("chr","",scPACds@anno$chr)
     sub_scPACds@anno$chr=gsub("chr","",sub_scPACds@anno$chr)
   }
    
      ##保留所有物种样本数据中主要染色体上的pA位点
      print("Filter PAC by main chrs")
      scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
      sub_scPACds=subsetPACds(sub_scPACds,chrs=chrs,minExprConds = 0,verbose =T)
    
    saveRDS(scPACds,file="PACds/raw.scPACds.rds")
    
    ##remove IP。
    scPACds =removePACdsIP(scPACds , bsgenome, returnBoth=FALSE, up=-140, dn=10, conA=6, sepA=NA,chrCheck = FALSE)
    
    #找出在100nt内具有polyA tails支持的pA位点
    scPACds@anno$PA_id=paste0(scPACds@anno$chr,":",scPACds@anno$strand,":",scPACds@anno$coord)
    sub_scPACds@anno$PA_id=paste0(sub_scPACds@anno$chr,":",sub_scPACds@anno$strand,":",sub_scPACds@anno$coord)
    
    #在原始的tail data 中去除IP位点
    info=intersect(scPACds@anno$PA_id,sub_scPACds@anno$PA_id)
    sub_scPACds@anno=sub_scPACds@anno[sub_scPACds@anno$PA_id %in% info,]
    sub_scPACds@counts=sub_scPACds@counts[row.names(sub_scPACds@anno),]
    
    scPACds@anno$pa.tails="no.tail"
    loc=match(sub_scPACds@anno$PA_id,scPACds@anno$PA_id)
    table(is.na(loc))
    
    scPACds@anno$pa.tails[loc]="tail"
    ##计算每个pa总的表达量
    scPACds@anno$count=Matrix::rowSums(scPACds@counts)
    
    #####注释PA位点
    scPACds<- annotatePAC(scPACds,gff)
    
    scPACds<- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
    
    
    data=subset(scPACds@anno,scPACds@anno$ftr=="3UTR")
    tail.data=subset(data,data$pa.tails=="tail")
    notail.data=subset(data,data$pa.tails=="no.tail")
    
    #sub_count=notail.data$count
    sub_count=tail.data$count
    
    #进行升序排列，并根据分位数（0.2，0.5，0.55，0.6）求分位值
    sub_count=sort(sub_count)
    count4=quantile(sub_count,probs=0.25,names=FALSE)
    count3=quantile(sub_count,probs=0.50,names=FALSE)
    count2=quantile(sub_count,probs=0.55,names=FALSE)
    count1=quantile(sub_count,probs=0.60,names=FALSE)
    
    #根据3UTR tail pa不同表达量的截断点进行分组
    scPACds@anno <-  scPACds@anno %>%
      mutate(count.group = case_when(
        count >count1 ~"count.rank1",
        count > count2 & count <= count1 ~ "count.rank2",
        count > count3 & count <= count2 ~ "count.rank3",
        count > count4 & count <= count3 ~ "count.rank4",
        count <=count4 ~"count.rank5"))
    
    scPACds@anno=subset(scPACds@anno,scPACds@anno$pa.tails=="tail" | scPACds@anno$count.group %in% c("count.rank1","count.rank2"))
    scPACds@counts=scPACds@counts[row.names(scPACds@anno),]
    
   
    print("Filter PAC by expression")
    scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
    summary(scPACds)
    
    #修改染色体名称：1-5，Mt，Pt--->Chr1-5,ChrM(线粒体染色体),ChrC
    #if(species=="arab"){
    #  scPACds@anno$chr=paste0("Chr",scPACds@anno$chr)}
    
    saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
    ##绘制3UTR区域的单核苷酸分布图
    
    fafiles = faFromPACds(scPACds, bsgenome
                          , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                          up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
    
   
    p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                          filepre = "")
    topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
    
    ######poly(A)信号
    #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
    # "ATTATA" "AACAAG" "AATAAG"
    #v=getVarGrams('mm')
    priority=c(1,2,rep(3, length(v)-2))
    
    scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
    scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                                 grams = v,
                                 priority = priority,
                                 from = -50, to = -1,
                                 label = "com", chrCheck = TRUE)
    
    pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
    
    pas$PAS[is.na(pas$PAS)]='NO PAS'
    #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
    pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
    
    message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    

    sink()
    cat(readLines("PACds/readPACs.log.txt"), sep="\n")
    
    #绘图
    n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
    n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
    
    n$PAC=n$nPAC/n2$nTot
    
    n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
    
    p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
      geom_bar(stat="identity",width = 0.7) +
      #scale_fill_brewer(palette="Set3" )+      
      scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
      ylab("PAS Fraction") + theme_classic()+
      geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
               
                vjust = -0.2, hjust = 0.5)+
      theme(legend.position = "none")
    topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
    
    rm(pas)
    
    pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
    p2=plotSeqLogo(pas)
    topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
    
    saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
    message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
  }
}

###process sierra PACs file---------------------------------------------------
#tool=tool.list[2]
toolfile="sierra"
tool="Sierra"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
 
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas

for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/",GSM)
  setwd(path)
  
 
  if (!dir.exists("./PACds")){
    dir.create("./PACds")
  } else {
    print("Output dir already exists!")
  }
  
  anno="peak_annotations.txt"
  count="peak.count"
  
  sink("PACds/readPACs.log.txt")
  #将pa读入PACdataset对象
  scPACds=unify.pa(anno,count,tool="Sierra")
  
   
    print("Filter PAC by main chrs")
    scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  
  saveRDS(scPACds,file="PACds/raw.scPACds.rds")
 
  print("Filter PAC by expression")
  scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
  
  ##remove IP:filtered false pA peaks based on the results of using `Sierra::AnnotatePeaksFromGTF` to analyze whether the pA peaks were from A- or T-rich regions.
  real.pa=subset(scPACds@anno,scPACds@anno$pA_stretch=="FALSE" & scPACds@anno$pT_stretch== "FALSE")
  scPACds=subsetPACds(scPACds, PAs = row.names(real.pa))
  print(paste0("Remove IP: subset real PAC(",nrow(scPACds@anno)," PACs)"))
  
  #####注释PA位点
  scPACds<- annotatePAC(scPACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(scPACds)
  
  saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(scPACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
  
  p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
  scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("PACds/readPACs.log.txt"), sep="\n")
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+      
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
             
              vjust = -0.2, hjust = 0.5)  
  topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
  message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
}
 }

###process scapture PACs file---------------------------------------------------
#tool=tool.list[3]
toolfile="scapture"
tool="SCAPTURE"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
 
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas

for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/",GSM)#when tools is scapture
  setwd(path)
  

  
  if (!dir.exists(paste0("./","PACds"))){
    dir.create(paste0("./","PACds"))
  } else {
    print("Output dir already exists!")
  }
  
  anno=paste0(GSM,".PASquant.KeepPAS.metadata")
  count=paste0(GSM,".PASquant.KeepCell.UMIs.tsv.gz")
  
  sink("PACds/readPACs.log.txt")
  #将pa读入PACdataset对象
  scPACds=unify.pa(anno,count,tool="SCAPTURE")
  
 
  print("Filter PAC by main chrs")
  scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(scPACds,file="PACds/raw.scPACds.rds")
  
  print("Filter PAC by expression")
  scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
  
  ##remove IP
  message("SCAPTURE has removed IP pA sites")
  
  #####注释PA位点
  scPACds<- annotatePAC(scPACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(scPACds)
  
  
  saveRDS(scPACds,file =paste0("./PACds/anno.scPACds.rds"))
  
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(scPACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
  #绘制3utr区域的PA位点的单核苷酸分布
  p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0("./PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
  scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("PACds/readPACs.log.txt"), sep="\n")
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+      
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
             
              vjust = -0.2, hjust = 0.5)+
    theme(legend.position = "none")   
  topptx(p, filename = paste0("./PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(scPACds3UTR,file = paste0("./PACds/scPACds3utr.rds"))
  message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
}
 }



###process Infernape PACs file---------------------------------------------------
#tool=tool.list[4]
toolfile="infernape"
tool="Infernape"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11

for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
  #不同数据下不同工具的输出文件路径
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas


for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/",GSM)
  setwd(path)
  

  if (!dir.exists("./PACds")){
    dir.create("./PACds")
  } else {
    print("Output dir already exists!")
  }
  
  anno="anno_filtered.csv"
  count="cnt_mat"
  
  sink("PACds/readPACs.log.txt")
  #将pa读入PACdataset对象
  
  if(species=="arab"){
    scPACds=unify.pa(anno,count=NULL,tool="Infernape")
    scPACds@anno$chr=gsub("Chr","",scPACds@anno$chr)
  }else if (species %in% c("human", "mouse")) {
    scPACds <- unify.pa(anno, count = count, tool = "Infernape")
  }
  
 
  print("Filter PAC by main chrs")
  scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(scPACds,file="PACds/raw.scPACds.rds")
  

  print("Filter PAC by expression")
  scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
  
  ##remove IP
  message("Infernape has removed IP pA sites")
  
  #####注释PA位点
  scPACds<- annotatePAC(scPACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(scPACds)
  
  
  saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(scPACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
  
  p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
  scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("PACds/readPACs.log.txt"), sep="\n")
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+     
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
             
              vjust = -0.2, hjust = 0.5)+
    theme(legend.position = "none")  
  topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
  message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
}
 }

###process SCAPE PACs file---------------------------------------------------
#tool=tool.list[5]
toolfile="scape"
tool="SCAPE"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)

  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas
  if(j==11){
    GSMlists=c("GSM3490690","GSM3490691")
  }


for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/",GSM)
  setwd(path)
  
 
  if (!dir.exists("./PACds")){
    dir.create("./PACds")
  } else {
    print("Output dir already exists!")
  }
  
  
  anno=paste0(file.path,"/",GSM,"/collapse_pa.tsv.gz")
  count=paste0(file.path,"/",GSM,"/pasite.csv.gz")
  
  sink("PACds/readPACs.log.txt")
  #将pa读入PACdataset对象
  scPACds=unify.pa(anno,count,tool="SCAPE")
  
 
  print("Filter PAC by main chrs")
  scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(scPACds,file="PACds/raw.scPACds.rds")
  
 
  print("Filter PAC by expression")
  scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
  
  ##remove IP
  message("SCAPE has removed IP pA sites")
  
  #####注释PA位点
  scPACds<- annotatePAC(scPACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(scPACds)
  
  
  saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
  ##绘制3UTR区域的单核苷酸分布图
 
  fafiles = faFromPACds(scPACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
 
  p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
  scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("PACds/readPACs.log.txt"), sep="\n")
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+     
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
              vjust = -0.2, hjust = 0.5)+
    theme(legend.position = "none")  
  topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
  message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
}
 }


###process scAPA PACs file---------------------------------------------------
#tool=tool.list[6]
toolfile="scapa"
tool="scAPA"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
  
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas

for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/",GSM)
  setwd(path)
  
 
  if (!dir.exists("./PACds")){
    dir.create("./PACds")
  } else {
    print("Output dir already exists!")
  }
  sink("PACds/readPACs.log.txt")
  
  if(species=="mouse" | species=="arab"){
    anno=paste0(file.path,"/",GSM,"/outs/out_peaks.txt")
    count=paste0(file.path,"/",GSM,"/outs/cell_Peaks.RDS")
    
    #将pa读入PACdataset对象
    scPACds=unify.pa(anno,count,tool="scAPA") 
  }
  
  if(species=="human"){
    message(paste0("read PACs from ", tool))
    anno=paste0(file.path,"/",GSM,"/outs/h38_peaks.bed")
    count=paste0(file.path,"/",GSM,"/outs/cell_Peaks.RDS")
    
    anno=read.table(anno,sep="\t",header = F)
    anno=anno[,c(1:4,6)]
    colnames(anno)=c("chr","start","end","geneid","strand")
    row.names(anno)=anno$geneid
    anno$coord <- 0
   
    anno[anno$strand == "+",]$coord <- anno[anno$strand == "+",]$end
   
    anno[anno$strand == "-",]$coord <- anno[anno$strand == "-",]$start
    
    anno=unique(anno)
    
    count=readRDS(count)
    ##exclude internal priming suspected peaks, peaks having a stretch of at least 8 consecutive As in the region between 10 nt to 140 nt. to the peak's end are filtered
    count=filter_IP_cell(x = count, int.priming.seq = "AAAAAAAA",
                      left = 10,
                      right = 140)
    
    
    #anno=anno[count@row.Data$GeneID,]
    
    #@cells.counts储存pa在单个细胞中的表达量
    count=count@cells.counts
    count=count[,-1]
    
    anno<-anno[intersect(row.names(count),row.names(anno)),]
    count<-count[intersect(row.names(count),row.names(anno)),]
    scPACds <- createPACdataset(counts=count, anno=anno)
    scPACds@colData$barcode<-row.names(scPACds@colData)
    #删掉group这一列
    scPACds@colData<- subset(scPACds@colData, select= barcode)
    print(paste0(length(row.names(scPACds@anno))," PACs","; ",length(colnames(scPACds@counts))," cells")) 
    
  }
  
 
  print("Filter PAC by main chrs")
  scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(scPACds,file="PACds/raw.scPACds.rds")
  
  ##remove IP:已在将pac数据读入pacdataset对象时去除
  message("scAPA has removed IP pA sites")
  
  
  print("Filter PAC by expression")
  scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
  
  #####注释PA位点
  scPACds<- annotatePAC(scPACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(scPACds)
  
  saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(scPACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
 
  p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
  scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("PACds/readPACs.log.txt"), sep="\n")
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+      
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
             
              vjust = -0.2, hjust = 0.5)+
    theme(legend.position = "none")  
  topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
  message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM))
  
}
 }

###process MAAPER PACs file---------------------------------------------------
#tool=tool.list[7]
toolfile="maaper"
tool="MAAPER"

celllists=c(2042,1827,8693,11175,5511,5316,510,3959,685,707,96,96,96,96,96,96,2471,4418,4476,4797,3119,5605,4406,1473,1643)
names(celllists)=c("merge","merge","merge","merge","MI_DAY7_GFP","SHAM_DAY7_GFP","GSM2906396","GSM2906399","GSM2333586","GSM2333587",
                   "GSM1524282","GSM1524283","GSM1524284","GSM1524285","GSM1524286","GSM1524287","3k","4k","5k","6k","GSM4712885","GSM4712907",
                   "GSM3490689","GSM3490690","GSM3490691")

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 1:3){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)

  file.path=paste0(GSEfile,"/",toolfile)
  
  #GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas

for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/",GSM)
  setwd(path)
  
  #细胞数
  if(j %in% 1:4){
    cell.num=celllists[j]
  }else{
    cell.num=celllists[names(celllists)==GSM]
  }
  
  
 
  if (!dir.exists("./PACds")){
    dir.create("./PACds")
  } else {
    print("Output dir already exists!")
  }
  
  if(j %in% 4:11){
    
    anno=paste0(file.path,"/",GSM,"/pas.txt")
    count=paste0(file.path,"/",GSM,"/gene.txt")
    
    sink("PACds/readPACs.log.txt") 
    
    
    #将pa读入PACdataset对象
    scPACds=unify.pa(anno,count,tool="MAAPER")
    
    
    print("Filter PAC by main chrs")
    scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
    
    saveRDS(scPACds,file="PACds/raw.scPACds.rds")
    
    
    min.tag=cell.num/200
    print("Filter PAC by expression")
    scPACds=subsetPACds(scPACds, totPACtag=min.tag,minExprConds=0,verbose =T)
    
  } else if(j %in% 1:3){
    pac=readRDS(paste0(path,"/paclists.rds"))
    
    sink("PACds/readPACs.log.txt") 
    scPACds=mergePACds_v0(pac,d=24)
    
    ##过滤低表达pA位点
    min.tag=cell.num/200
    filter_edgeR <- edgeR::filterByExpr(scPACds@counts,
                                        design = NULL,
                                        group = scPACds@colData$celltype,
                                        lib.size = NULL,
                                        min.count = 0,
                                        min.total.count =min.tag,
                                        large.n = 20,
                                        min.prop = 0.5) 
    table(filter_edgeR)
    scPACds@counts= scPACds@counts[filter_edgeR,]
    scPACds@anno= scPACds@anno[filter_edgeR,]
    
  }
  
  
  ##remove IP
  message("MAAPER has removed IP pA sites")
  
  #####注释PA位点
  scPACds<- annotatePAC(scPACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(scPACds)
  
  saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(scPACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
  
  p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
  scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
 cat(readLines("PACds/readPACs.log.txt"), sep="\n")
 
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+     
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
              vjust = -0.2, hjust = 0.5) +
    theme(legend.position = "none") 
  topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
  message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
}
 }

###process polyapipe PACs file---------------------------------------------------
#tool=tool.list[8]
toolfile="polyapipe"
tool="polyApipe"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
 
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas
  
  for(i in 1: length(GSMlists)){
    message("start") 
    GSM=GSMlists[i]
    path=paste0(file.path,"/",GSM)
    setwd(path)
    
   
    if (!dir.exists("./PACds")){
      dir.create("./PACds")
    } else {
      print("Output dir already exists!")
    }
    
    anno=paste0(GSM,"_polyA_peaks.gff")
    count=paste0(GSM,"_counts.tab.gz")
    
    sink("PACds/readPACs.log.txt")
    #将pa读入PACdataset对象
    scPACds=unify.pa(anno,count,tool="polyApipe")
    
   
    print("Filter PAC by main chrs")
    scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
    
    
    saveRDS(scPACds,file="PACds/raw.scPACds.rds")
    
   
    print("Filter PAC by expression")
    scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
    
    ##remove IP
   
    real.pa=subset(scPACds@anno,scPACds@anno$misprime=="False")
    scPACds=subsetPACds(scPACds, PAs = row.names(real.pa))
    print(paste0("Remove IP: subset real PAC(",nrow(scPACds@anno)," PACs)"))
    
    #####注释PA位点
    scPACds<- annotatePAC(scPACds,gff)
    
    ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
    scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
    summary(scPACds)
    
    saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
    ##绘制3UTR区域的单核苷酸分布图
    
    fafiles = faFromPACds(scPACds, bsgenome
                          , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                          up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
    
    
    p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                          filepre = "")
    topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
    
    ######poly(A)信号的计算
    #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
    # "ATTATA" "AACAAG" "AATAAG"
    #v=getVarGrams('mm')
    priority=c(1,2,rep(3, length(v)-2))
    
    scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
    scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                                 grams = v,
                                 priority = priority,
                                 from = -50, to = -1,
                                 label = "com", chrCheck = TRUE)
    
    pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
    
    pas$PAS[is.na(pas$PAS)]='NO PAS'
    #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
    pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
    
    message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    
    print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    
    sink()
    cat(readLines("PACds/readPACs.log.txt"), sep="\n")
    
    #绘图
    n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
    n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
    
    n$PAC=n$nPAC/n2$nTot
    
    n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
    
    p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
      geom_bar(stat="identity",width = 0.7) +
      #scale_fill_brewer(palette="Set3" )+      
      scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
      ylab("PAS Fraction") + theme_classic()+
      geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
                vjust = -0.2, hjust = 0.5)  
    topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
    
    rm(pas)
    
    pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
    p2=plotSeqLogo(pas)
    topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
    
    saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
    message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
  }
}

###process scraps PACs file---------------------------------------------------
#tool=tool.list[9]
toolfile="scraps"
tool="scraps"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 11:11){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
 
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas
  
  for(i in 1: length(GSMlists)){
    message("start") 
    GSM=GSMlists[i]
    path=paste0(file.path,"/",GSM)
    setwd(path)
    
    
    if (!dir.exists("./PACds")){
      dir.create("./PACds")
    } else {
      print("Output dir already exists!")
    }
    
    if(j %in% 1:4){
      pafile=paste0(GSE,"_R2_counts.tsv.gz")#when j=1:4
    }
    if(j %in% c(5,7:11)){
      pafile=paste0(GSM,"_R2_counts.tsv.gz")#when j=5:11
    }
    
    if(j== 6){
      pafile=paste0(GSM,"_paired_counts.tsv.gz")#when j=6,mouse_bone sample
    }
    
    sink("PACds/readPACs.log.txt")
    #将pa读入PACdataset对象
    scPACds=unify.pa(pafile=pafile,tool="scraps")
    
   
    print("Filter PAC by main chrs")
    scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
    
    
    saveRDS(scPACds,file="PACds/raw.scPACds.rds")
    
    print("Filter PAC by expression")
    scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
    
    ##remove IP:使用了polyadb中的pA位点做参考，可不再过滤
    message("scraps has removed IP pA sites")
    
    
    #####注释PA位点
    scPACds<- annotatePAC(scPACds,gff)
    
    ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
    scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
    summary(scPACds)
    
    saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
    ##绘制3UTR区域的单核苷酸分布图
    
    fafiles = faFromPACds(scPACds, bsgenome
                          , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                          up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
    
    
    p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                          filepre = "")
    topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
    
    ######poly(A)信号的计算
    #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
    # "ATTATA" "AACAAG" "AATAAG"
    #v=getVarGrams('mm')
    priority=c(1,2,rep(3, length(v)-2))
    
    scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
    scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                                 grams = v,
                                 priority = priority,
                                 from = -50, to = -1,
                                 label = "com", chrCheck = TRUE)
    
    pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
    
    pas$PAS[is.na(pas$PAS)]='NO PAS'
    #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
    pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
    
    message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    
    print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    
    sink()
    cat(readLines("PACds/readPACs.log.txt"), sep="\n")
    
    #绘图
    n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
    n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
    
    n$PAC=n$nPAC/n2$nTot
    
    n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
    
    p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
      geom_bar(stat="identity",width = 0.7) +
      #scale_fill_brewer(palette="Set3" )+      
      scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+
      ylab("PAS Fraction") + theme_classic()+
      geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
                vjust = -0.2, hjust = 0.5)  
    topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
    
    rm(pas)
    
    pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
    p2=plotSeqLogo(pas)
    topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
    
    saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
    message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
  }
}


###process scUTRquant pac files----------------------------------------------------


#tool=tool.list[10]
toolfile="scutrquant"
tool="scUTRquant"

#when species is mouse,j=1:8;specise is human,j=9:10;specise is arab,j=11
for(j in 9:10){
  GSE=GSElists[j]
  GSEfile=paste0("G:/scAPA/PAC/",GSE)
  
  file.path=paste0(GSEfile,"/",toolfile)
  
  GSMlists=dir(file.path)#whenGSE=GSElists[5:11]
  #GSMlists="merge" #whenGSE=GSElists[1:4]
  #GSMlists=datas
  
  for(i in 1: length(GSMlists)){
    message("start") 
    GSM=GSMlists[i]
    path=paste0(file.path,"/",GSM)
    setwd(path)
    
    
    if (!dir.exists("./PACds")){
      dir.create("./PACds")
    } else {
      print("Output dir already exists!")
    }
    
      pafile=paste0(GSM,".txs.rds")#when j=5:11
    
    
    sink("PACds/readPACs.log.txt")
    scPACds=unify.pa(pafile=pafile,tool="scUTRquant")
    
    
    print("Filter PAC by main chrs")
    scPACds=subsetPACds(scPACds,chrs=chrs,minExprConds = 0,verbose =T)
    
    
    saveRDS(scPACds,file="PACds/raw.scPACds.rds")
    
    
    print("Filter PAC by expression")
    scPACds=subsetPACds(scPACds, totPACtag=1,minExprConds=ncol(scPACds@counts)/200,verbose =T)
    
    ##remove IP:使用了polyadb中的pA位点做参考，可不再过滤
    message("scUTRquant has removed IP pA sites")
    
    
    #####注释PA位点
    scPACds<- annotatePAC(scPACds,gff)
    
    ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
    scPACds <- ext3UTRPACds(scPACds,ext3UTRlen = ext.len,extFtr='3UTR')
    summary(scPACds)
    
    saveRDS(scPACds,file = "PACds/anno.scPACds.rds")
    ##绘制3UTR区域的单核苷酸分布图
    
    fafiles = faFromPACds(scPACds, bsgenome
                          , what = "updn", fapre = paste0(tool,"_",genome,".200nt"),
                          up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
    
  
    p=plotATCGforFAfile.2(paste0(tool,"_",genome,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                          filepre = "")
    topptx(p, filename =paste0("PACds/",tool,"_",genome,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
    
    ######poly(A)信号的计算
    #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
    # "ATTATA" "AACAAG" "AATAAG"
    #v=getVarGrams('mm')
    priority=c(1,2,rep(3, length(v)-2))
    
    scPACds3UTR<-scPACds[scPACds@anno$ftr=='3UTR']
    scPACds3UTR <- annotateByPAS(pacds = scPACds3UTR, bsgenome,
                                 grams = v,
                                 priority = priority,
                                 from = -50, to = -1,
                                 label = "com", chrCheck = TRUE)
    
    pas=as.data.frame(cbind(region=scPACds3UTR@anno$ftr, PAS=scPACds3UTR@anno$com_gram))
    
    pas$PAS[is.na(pas$PAS)]='NO PAS'
    #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
    pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
    
    message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    
    print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
    
    sink()
    cat(readLines("PACds/readPACs.log.txt"), sep="\n")
    
    #绘图
    n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
    n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
    
    n$PAC=n$nPAC/n2$nTot
    
    n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
    
    p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
      geom_bar(stat="identity",width = 0.7) +
      #scale_fill_brewer(palette="Set3" )+      
      scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+
      ylab("PAS Fraction") + theme_classic()+
      geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'),
                vjust = -0.2, hjust = 0.5)  
    topptx(p, filename = paste0("PACds/",tool,"_3UTR.PAS.pptx"), width = 6, height = 4)
    
    rm(pas)
    
    pas=scPACds3UTR@anno$com_gram[!is.na(scPACds3UTR@anno$com_gram)]
    p2=plotSeqLogo(pas)
    topptx(p2, filename = paste0("PACds/",tool,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
    
    saveRDS(scPACds3UTR,file = "PACds/scPACds3utr.rds")
    message(paste0(toolfile,": ",i,"/",length(GSMlists),GSM)) 
  }
}


###process bulk pac data-------------------------------------------------
#######(1)polyAseqTrap################----------------------------------------------------------

tool="polyAseqTrap"


GSE="mouse_esc"
file.path=paste0("G:/scAPA/bulk/",GSE)
setwd(file.path)


datas=c("SRR1033836")
GSMlists=datas
min.tag=2

paclists=list()

for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/PACds")
  #setwd(path)
  
 
  if (!dir.exists(path)){
    dir.create(path)
  } else {
    print("Output dir already exists!")
  }
  
  output.path=paste0(path,"/",GSM)
  if (!dir.exists(output.path)){
    dir.create(output.path)}
  
  setwd(output.path)
  
  pafile=paste0(file.path,"/",GSM,"_PAtable_renew.Rdata")
  
 
  sink("readPACs.log.txt") 
  
  PACds=unify.pa(pafile=pafile ,tool="polyAseqTrap")
  
  
  print("Filter PAC by main chrs")
  PACds=subsetPACds(PACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(PACds,file="raw.PACds.rds")
  
  
  print("Filter PAC by expression")
  PACds=subsetPACds(PACds, totPACtag=min.tag,minExprConds=0,verbose =T)
  
  #####注释PA位点
  PACds<- annotatePAC(PACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  PACds <- ext3UTRPACds(PACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(PACds)
  
  saveRDS(PACds,file="anno.PACds.rds")
  paclists[GSM]=PACds
  
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(PACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",GSM,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
 
  p=plotATCGforFAfile.2(paste0(tool,"_",GSM,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0(tool,"_",GSM,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  PACds3UTR<-PACds[PACds@anno$ftr=='3UTR']
  PACds3UTR <- annotateByPAS(pacds = PACds3UTR, bsgenome,
                               grams = v,
                               priority = priority,
                               from = -50, to = -1,
                               label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=PACds3UTR@anno$ftr, PAS=PACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("readPACs.log.txt"), sep="\n")
  
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+      
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
              vjust = -0.2, hjust = 0.5) +
    theme(legend.position = "none") 
  topptx(p, filename = paste0(tool,"_",GSM,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=PACds3UTR@anno$com_gram[!is.na(PACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0(tool,"_",GSM,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(PACds3UTR,file = "PACds3utr.rds")
  message(paste0(tool,": ",i,"/",length(GSMlists),GSM)) 
}

saveRDS(paclists,file=paste0(file.path,"/PACds/paclists.rds"))

#######(2)RAW.PAC##############--------------------------------
tool="RWA.PAC"

GSE="mouse_esc_mef"
file.path=paste0("G:/scAPA/bulk/",GSE,"/RAW.PAC")
#file.path=paste0("G:/scAPA/bulk/",GSE,"/RAW.PAC/polyasite")
setwd(file.path)


datas=c("mouse_esc_mef")
GSMlists=datas
min.tag=2

for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/PACds")
  #setwd(path)
  

  if (!dir.exists(path)){
    dir.create(path)
  } else {
    print("Output dir already exists!")
  }
  
  setwd(path)
  
  sink("readPACs.log.txt") 
  
  
  PACds=readRDS(paste0(file.path,"/paclists.rds"))
 
  PACds=mergePACds(PACds,d=24)
  
  print(paste0(length(row.names(PACds@anno))," PACs"))
  
  print("Filter PAC by main chrs")
  PACds=subsetPACds(PACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(PACds,file="raw.PACds.rds")
  
  
  print("Filter PAC by expression")
  PACds=subsetPACds(PACds, totPACtag=min.tag,minExprConds=0,verbose =T)
  
  #####注释PA位点
  PACds<- annotatePAC(PACds,gff)
  
  ####向3utr向下游延申，将位于其中的PA位点归为3utrPA位点
  PACds <- ext3UTRPACds(PACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(PACds)
  
  saveRDS(PACds,file="anno.PACds.rds")
 
  
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(PACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",GSM,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
  
  p=plotATCGforFAfile.2(paste0(tool,"_",GSE,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0(tool,"_",GSE,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  PACds3UTR<-PACds[PACds@anno$ftr=='3UTR']
  PACds3UTR <- annotateByPAS(pacds = PACds3UTR, bsgenome,
                             grams = v,
                             priority = priority,
                             from = -50, to = -1,
                             label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=PACds3UTR@anno$ftr, PAS=PACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("readPACs.log.txt"), sep="\n")
  
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+     
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
  
              vjust = -0.2, hjust = 0.5) +
    theme(legend.position = "none") 
  topptx(p, filename = paste0(tool,"_",GSE,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=PACds3UTR@anno$com_gram[!is.na(PACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0(tool,"_",GSE,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(PACds3UTR,file = "PACds3utr.rds")
  message(paste0(tool,": ",i,"/",length(GSMlists),GSM)) 
}


####(3)QAPA#######-------------------------------------------------------------

tool="QAPA"


GSE="mouse_tip"
file.path=paste0("G:/scAPA/bulk/",GSE)
setwd(file.path)


datas=c("mouse_tip")
GSMlists=datas
##mouse_tip:共48个样本
min.tag=24



for(i in 1: length(GSMlists)){
  message("start") 
  GSM=GSMlists[i]
  path=paste0(file.path,"/PACds")
  #setwd(path)
  
  
  if (!dir.exists(path)){
    dir.create(path)
  } else {
    print("Output dir already exists!")
  }
  
  #output.path=paste0(path,"/",GSM)
  #if (!dir.exists(output.path)){
  #  dir.create(output.path)}
  
  setwd(path)
  
  pafile=paste0(file.path,"/","qapa.results.txt")
  
  #记录控制台的输出内容
  sink("readPACs.log.txt") 
  
  
  #将pa读入PACdataset对象
  PACds=unify.pa(pafile=pafile ,tool="QAPA")
  
  
  print("Filter PAC by main chrs")
  PACds=subsetPACds(PACds,chrs=chrs,minExprConds = 0,verbose =T)
  
  saveRDS(PACds,file="raw.PACds.rds")
  
  
  PACds=mergePACds(PACds,d=24)
  

  print("Filter PAC by expression")
  PACds=subsetPACds(PACds, totPACtag=min.tag,minExprConds=0,verbose =T)
  
  #####注释PA位点
  PACds<- annotatePAC(PACds,gff)
  
  PACds <- ext3UTRPACds(PACds,ext3UTRlen = ext.len,extFtr='3UTR')
  summary(PACds)
  
  saveRDS(PACds,file="anno.PACds.rds")
  
  
  ##绘制3UTR区域的单核苷酸分布图
  
  fafiles = faFromPACds(PACds, bsgenome
                        , what = "updn", fapre = paste0(tool,"_",GSM,".200nt"),
                        up = -100, dn = 100, byGrp = "ftr", chrCheck = FALSE)
  
 
  p=plotATCGforFAfile.2(paste0(tool,"_",GSM,".200nt.3UTR.fa"), ofreq = FALSE, opdf = FALSE, refPos = 101,
                        filepre = "")
  topptx(p, filename =paste0(tool,"_",GSM,".200nt.3UTR.pptx"), width = 6, height = 4,units = "in")
  
  ######poly(A)信号的计算
  #人和鼠周围常见的poly(A)信号: "AATAAA" "ATTAAA" "TATAAA" "AGTAAA" "AATACA" "CATAAA" "AATATA" "GATAAA" "AATGAA" "AATAAT" "AAGAAA" "ACTAAA" "AATAGA" "ATTACA" "AACAAA"
  # "ATTATA" "AACAAG" "AATAAG"
  #v=getVarGrams('mm')
  priority=c(1,2,rep(3, length(v)-2))
  
  PACds3UTR<-PACds[PACds@anno$ftr=='3UTR']
  PACds3UTR <- annotateByPAS(pacds = PACds3UTR, bsgenome,
                             grams = v,
                             priority = priority,
                             from = -50, to = -1,
                             label = "com", chrCheck = TRUE)
  
  pas=as.data.frame(cbind(region=PACds3UTR@anno$ftr, PAS=PACds3UTR@anno$com_gram))
  
  pas$PAS[is.na(pas$PAS)]='NO PAS'
  #pas$PAS[!pas$PAS %in% c("AATAAA","ATTAAA","AAGAAA","NOPAS")]='other PAS'
  pas$PAS[!pas$PAS %in% c("AATAAA","NO PAS")]='1nt variants'
  
  message(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  print(paste0("PAS.number: ", length(pas$PAS[!pas$PAS=="NO PAS"])))
  
  sink()
  cat(readLines("readPACs.log.txt"), sep="\n")
  
  
  #绘图
  n=pas %>% dplyr::group_by(region, PAS)  %>% dplyr::summarise(nPAC=n())
  n2=pas %>% dplyr::group_by(region)  %>% dplyr::summarise(nTot=n())
  
  n$PAC=n$nPAC/n2$nTot
  
  n$PAS=factor(n$PAS, levels=rev(c("NO PAS","1nt variants","AATAAA")))
  
  p=ggplot(data=n, aes(x=PAS, y=PAC, fill=PAS)) +
    geom_bar(stat="identity",width = 0.7) +
    #scale_fill_brewer(palette="Set3" )+     
    scale_fill_manual(values= c("#299D8F","#E9C46A","#D87659"))+ 
    ylab("PAS Fraction") + theme_classic()+
    geom_text(aes(label=round(PAC,2)), position = position_dodge2(width = 0.9, preserve = 'single'), 
              vjust = -0.2, hjust = 0.5) +
    theme(legend.position = "none") 
  topptx(p, filename = paste0(tool,"_",GSM,"_3UTR.PAS.pptx"), width = 6, height = 4)
  
  rm(pas)
  
  pas=PACds3UTR@anno$com_gram[!is.na(PACds3UTR@anno$com_gram)]
  p2=plotSeqLogo(pas)
  topptx(p2, filename = paste0(tool,"_",GSM,"_3UTR.PAS_logo.pptx"), width = 6, height = 4)
  
  saveRDS(PACds3UTR,file = "PACds3utr.rds")
  message(paste0(tool,": ",i,"/",length(GSMlists),GSM)) 
}




