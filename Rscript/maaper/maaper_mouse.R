#maaper requires three input files:

#    The GTF file of the reference genome;
#    The BAM files of the 3’ sequencing data (nearSite reads). The BAM file should be sorted and the index BAI file should be present in the same directory as the BAM file;
#    The PAS annotation file whose version matches the reference genome. We have prepared PolyA_DB annotation files for MAAPER, and they can be downloaded from this page.

#The final output of mapper are two text files named “gene.txt” and “pas.txt”, which contain the predicted PASs and APA results.

rm(list = ls())

library(MAAPER)

##Record the memory used and the time
start.time=Sys.time()
start_mem <- as.numeric(system(paste("ps -p", Sys.getpid(), "-o rss="), intern = TRUE)) / 1024


pas_annotation = readRDS("/mnt/user/BXY/scAPA/scAPA_tools/supplementary_data/MAAPER/PolyA_DB/mouse.PAS.mm10.rds")
gtf = "/mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf"

## input path----
input_path='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/mouse_sperm'

## output dir----
out_path='/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/mouse_sperm/maaper'
GSMlists=c("pas1_gn5000_bn2000_rep1","pas1_gn5000_bn2000_rep2","pas1_gn5000_bn2000_rep3")

##read2 length
read_lens=c(98,98,98)

for(i in 1:length(GSMlists)){
GSM=GSMlists[i]
read_len=read_lens[i]
##NOPAIRED TEST--------------------------------------
###########find peaks from all cells, not cluster
# # bam file of condition 1 /condition 2 (could be a vector if there are multiple samples)

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
