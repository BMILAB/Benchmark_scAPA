
rm(list = ls())

library(Infernape)

##Record runtime and memory usage
start.time=Sys.time()
start_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024

####infernape:find and count peaks

###input file path------
file.path="/mnt/T20T-6d12/BXY/scAPA/sim_data/SRA/human_GSM4712885"

GSMlists=c("pas1_gn5000_bn3000_rep1","pas1_gn5000_bn3000_rep2","pas1_gn5000_bn3000_rep3")

for (i in 1: length(GSMlists)){

GSM=GSMlists[i]

input.path=file.path
bam=paste0(input.path,"/",GSM,".bam")

##Read the cellbarcode and modify its format
##For the cellbarcode read by Infernape, the first column is the row number, the second column is labeled “x”, and the data is the barcode
whitelist.file=read.delim2(paste0(input.path,"/",GSM,".barcode_list.txt"),header=FALSE)$V1

utils::write.csv(whitelist.file,file=paste0(input.path,"/",GSM,"_infernape.bc.csv"),row.names = TRUE,quote=FALSE)

barcode=paste0(input.path,"/",GSM,"_infernape.bc.csv")

###output file path-----
output.path=paste0("/mnt/T20T-6d12/BXY/scAPA/sim_data/PAC/human_GSM4712885/infernape","/",GSM)

if(!dir.exists(output.path)){
    dir.create(output.path)
    }else{
    print("Dir already exists")
    }


##genome file------
##huamn:BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
##mouse:BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
genome = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38


##gene annotation file---------------------------
genome.ref="/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/GRch38fa+gtf/genecode/v34/GRch38_infernape_ref.csv"

gene.start=1
gene.end=2912496

pas.ref="/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/scAPA_tools/Infernape/Infernape_reference/Human/anno.polya_db.hg38.csv"

Infernape_cnt(genome.ref = genome.ref, #gtf file
              bam = bam, 
              batch.start = gene.start, #ID of the first gene in this batch
              batch.end = gene.end, #ID of the last gene in this batch
              ncores = 12, 
              d = 31,   #A parameter for smoothing
              h = 5,    #A parameter for smoothing
              d.cut = 50, #A parameter for merging close peaks
              hr = 160,   #The estimated half range of peak mode - PAS interval. Used for annotating peak modes
              min.mode.prop = 0.05, #Minimal relative threshold (wrt highest peak in 3'UTR region) to filter out small peaks
              min.mode.cutoff = 5,  #Minimal absolute coverage threshold to filter out small peaks
              output.path = output.path,
              pas.reference.file = pas.ref, #Known PAS table
              genome = genome,                                 #Full genome sequence. Required to be a S4 object.
              pas.search.cut.1 = 0, #Distance (nt) from the peak mode to upstream boundary of the searching interval for PAS
              pas.search.cut.2 = 300, #Distance (nt) from the peak mode to downstream boundary of the searching interval for PAS
              polystretch_length = 13, #Length of consecutive A sequence
              max_mismatch = 1,        #Maximal tolerance of mismatch
              motif.search.cut = 300,  #Window width for searching specified motifs
              invert_strand = FALSE,
              q = c(110, 200),   #Vector of length 2 which defines the PAS searching window
              whitelist.file = barcode,
              start.cid = NULL,  #The first peak cluster ID to analyze
              end.cid = NULL     #The last peak cluster ID to analyze
)
}

#Record runtime and memory usage
end.time=Sys.time()

end_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024
peak_mem <- system(
  sprintf("ps -o rss= -p %d --sort=-rss | head -1", Sys.getpid()),
  intern = TRUE) |> as.numeric() / 1024

cat("\n===== Infernape 性能分析报告 =====\n")
print(end.time-start.time)  

cat(sprintf("内存使用:\n- 初始: %.1f MB\n- 结束: %.1f MB\n- 峰值: %.1f MB\n",
            start_mem, end_mem, peak_mem))

cat("========================\n")

