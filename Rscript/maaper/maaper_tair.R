
rm(list = ls())

library(MAAPER)

##Record the memory used and the time
start.time=Sys.time()
start_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024


pas_annotation = readRDS("/mnt/user/BXY/ref_file/refPAC/tair10/maaper_tair10_pas.rds")
gtf = "/mnt/user/BXY/ref_file/TAIR/Arabidopsis_thaliana.TAIR10.49.gtf"

## input path----
input_path='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/tair_GSM3490690'

## output dir----
out_path='/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/tair_GSM3490690/maaper'
GSMlists=c("pas1_gn1000_bn1400_rep1","pas1_gn1000_bn1400_rep2","pas1_gn1000_bn1400_rep3")

##read2 length
read_lens=c(115,115,115)

for(i in 1:length(GSMlists)){
GSM=GSMlists[i]
read_len=read_lens[i]
##NOPAIRED TEST--------------------------------------
###########find peaks from all cells, not cluster
# bam file of condition 1 /condition 2 (could be a vector if there are multiple samples)

bam_c1 = paste0(input_path,"/",GSM, ".bam")


# output directory( not need create folder in advance)
output_dir = paste0(out_path,"/" ,GSM)
if (!dir.exists(output_dir))
        dir.create(output_dir)

maaper(gtf, # full path of the GTF file
       pas_annotation, # PAS annotation
       output_dir = output_dir, 
       bam_c1, bam_c1, # full path of the BAM files
       read_len = read_len, # read length
       ncores = 12,  # number of cores used for parallel computation 
       paired = FALSE,
       bed = TRUE# obtain bedGraph files for visualization with UCSC or IGV genome browser
      )
}   
   
# Record the memory used and the time
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
