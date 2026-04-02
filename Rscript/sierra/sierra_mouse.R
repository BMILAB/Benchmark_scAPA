####********Sierra*************######
#regtools junctions extract -s RF filterCB.bam -o data_junctions.bed
rm(list = ls())


library(BSgenome.Mmusculus.UCSC.mm10)
library(Sierra)

##Record the memory used and the time
start.time=Sys.time()
start_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024


## input path----
input_path='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/mouse_sperm'

## output dir----
out_path='/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/mouse_sperm/sierra'
GSMlists=c("pas1_gn5000_bn2000_rep1","pas1_gn5000_bn2000_rep2","pas1_gn5000_bn2000_rep3")

for(i in 1:length(GSMlists)){
GSM=GSMlists[i]

whitelist.file<- paste0(input_path,"/",GSM,'.barcode_list.txt') 
#whitelist.file<- paste0(input_path,'/barcode.txt') 
reference.file <-'/mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf' 
junctions.file <-paste0(input_path,"/",GSM,'_junctions.bed')
bamfile <- paste0(input_path,"/",GSM, ".bam")


output_dir <- paste0(out_path,"/",GSM)
# print(output_dir)
if (!dir.exists(output_dir)){
    dir.create(output_dir)
} else {
    print("Dir already exists!")
} 

setwd(output_dir)

##1.peak calling----------------------------------------------------------------------------

chrs=paste0("chr", c(as.character(1:19), 'X' ,'Y'))
FindPeaks(output.file ="peaks.txt",    # output filename
          chr.names = chrs ,
          filter.chr =TRUE,
          gtf.file = reference.file,              # gene model as a GTF file
          bamfile =bamfile,                  # BAM alignment filename.
          junctions.file = junctions.file,         #bed file     
          ncores = 12) 
      
peak<-read.table("peaks.txt", header = T, sep = "\t", stringsAsFactors = FALSE,fill=TRUE)


peak<-peak[peak$polyA_ID!="",]
all.peaks <- peak$polyA_ID
strand = data.frame(name=rownames(peak),s=sub(".*:.*:.*-.*:(.*)", "\\1", all.peaks))
strand<-strand[strand$s %in% c(1,-1),]
peak<-peak[as.character(strand$name),]


output.file<-"extract_peaks.txt"
write.table(peak, file = output.file, 
          quote = FALSE, row.names = FALSE,sep="\t")
          


##2.count peak--------------------------------------------------------------------------------------------------------------------------------------
#Generates a UMI count matrix where rows are the peaks and columns are the cells.
#Counts cells that are identified through a provided 'white list' of cell barcodes.

output.file<-"extract_peaks.txt"

count.dirs<-"peak.count" #not need creat folder in advance
CountPeaks(peak.sites.file = output.file, 
           gtf.file = reference.file,
           bamfile = bamfile, 
           whitelist.file = whitelist.file,
           output.dir = count.dirs, 
           countUMI = TRUE, 
           ncores = 12)


##3.annotate peak-------------------------------------------------------------------------------------------------------------------------------
## The annotation step will identify if peaks overlap exons, introns or UTRs. 
##AnnotatePeaksFromGTF takes as input the peak coordinates file, GTF file and genome file(optional )in addition to an output file name. 
#pA_motif_max_position:Any AAUAAA after this position are not considered (default 50nt)
#AAA_motif_min_position:Any polyA/polyT stretches before this postion are not considered (default 10)
#BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

bsgenome=BSgenome.Mmusculus.UCSC.mm10
AnnotatePeaksFromGTF(peak.sites.file = output.file,
                     genome=bsgenome,
                     pA_motif_max_position=50,
                     AAA_motif_min_position=10,
                     gtf.file = reference.file,
                     output.file = "peak_annotations.txt")

}

##Record the memory used and the time
end.time=Sys.time()

end_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024
peak_mem <- system(
  sprintf("ps -o rss= -p %d --sort=-rss | head -1", Sys.getpid()),
  intern = TRUE) |> as.numeric() / 1024

cat("\n===== 性能分析报告 =====\n")
print(end.time-start.time)  

cat(sprintf("内存使用:\n- 初始: %.1f MB\n- 结束: %.1f MB\n- 峰值: %.1f MB\n",
            start_mem, end_mem, peak_mem))

cat("========================\n")

 
