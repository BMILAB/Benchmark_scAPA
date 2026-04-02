rm(list = ls())
library(scAPAtrap)

##Record the memory used and the time
start.time=Sys.time()

start_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024

library(scAPAtrap)

start.time=Sys.time()
## full path of tools-----
tools=list(samtools='/home/test/miniconda3/envs/R4.3/bin/samtools',
           umitools='/home/test/miniconda3/envs/R4.3/bin/umi_tools',
           featureCounts="/home/test/miniconda3/envs/R4.3/bin/featureCounts",
           star='/home/test/miniconda3/envs/R4.3/bin/STAR')
           
## input path----
dir0='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/human_GSM4712885'

## output dir----
out_path="/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/human_GSM4712885/scapatrap"
GSMlists=c("pas1_gn5000_bn3000_rep1","pas1_gn5000_bn3000_rep2","pas1_gn5000_bn3000_rep3")

#read2 length
read_lens=c(91,91,91)

for(i in 1:length(GSMlists)){
GSM=GSMlists[i]

#inputBam=paste0(dir0,"/",GSM, "/filterCB.bam")
inputBam=paste0(dir0,"/",GSM, ".bam")

## log file to LOG all information (time, command, output file names..)
logf=gsub('.bam', '.APA.tails.onestep.log', inputBam, fixed=TRUE)


#output path will be under inputBam's dir if only dirname is provided(only provided name,not created); otherwise use the full path)
output_dir <- paste0(out_path,"/",GSM)
#output_dir <- paste0(out_path)

# print(output_dir)
if (!dir.exists(output_dir)){
   dir.create(output_dir)
} else {
    print("Dir already exists!")
} 

outputDir=paste0(output_dir,"/APA.tails")
#outputDir=paste0(output_dir)


##set parameters------
trap.params=setTrapParams()

#trap.params$barcode <- read.delim2(paste0(dir0, "/",GSM,'/barcode.txt'), header = F)$V1
trap.params$barcode <- read.delim2(paste0(dir0, "/",GSM,'.barcode_list.txt'), header = F)$V1

#delete "-1" of barcode
trap.params$barcode <-gsub("-[0-9]","",trap.params$barcode)
## set tail search way("genome","peaks","no")
trap.params$tails.search='genome'
#chr.name
trap.params$chrs=paste0("chr",c(as.character(1:22),'X','Y'))
#trap.params$chrs=c(as.character(1:5))
#bam contain CB,UB?
trap.params$TenX=TRUE

trap.params$maxwidth=1000
#read2 length
trap.params$readlength=read_lens[i]

trap.params$cov.cutoff=2

trap.params$min.cells=2

trap.params$min.count=2
trap.params$thread=12

## Run scAPAtrap
scAPAtrap(tools=tools,
          trap.params=trap.params,
          inputBam=inputBam,
          outputDir=outputDir,
          logf=logf)
}          
##check pas result
#load("/mnt/user/BXY/scAPA/PAC/scAPAtrap/GSM2803334/APA.notails/scAPAtrapData.rda")
#scAPAtrapData$peaks.count[,1:8]
#head(scAPAtrapData$peaks.meta)

#Record the memory used and the time
end.time=Sys.time()

end_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024
peak_mem <- system(
  sprintf("ps -o rss= -p %d --sort=-rss | head -1", Sys.getpid()),
  intern = TRUE) |> as.numeric() / 1024

cat("\n===== scAPAtrap性能分析报告 =====\n")
print(end.time-start.time)  

cat(sprintf("内存使用:\n- 初始: %.1f MB\n- 结束: %.1f MB\n- 峰值: %.1f MB\n",
            start_mem, end_mem, peak_mem))

cat("========================\n")


















