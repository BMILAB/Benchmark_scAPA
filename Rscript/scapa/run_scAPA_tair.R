##scAPA:Defining 3'UTR peaks and Quantifying the usage of each peak
# R code
require(package = "dplyr", warn.conflicts = F)
require(package = "tidyr", warn.conflicts = F)
require(package = "ggplot2", warn.conflicts = F)
require(package = "EnvStats", warn.conflicts = F)
require(package = "parallel", warn.conflicts = F)
require(package = "mclust", warn.conflicts = F)
require("Rsubread")
require("scAPA")
source("/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/code/scAPA/read_down.seq2.R")
source("/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/code/scAPA/intersect_peaks2.R")

##Record runtime
start.time=Sys.time()

####fasta files----------------------------------------------
char.length.path="/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/chrNameLength.txt"
fasta.path="/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

##When sample is Arabispis,need 3UTR bed file
utr="/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/scapa_utr3.bed"

###read cell annotation file----------------------------------
samplelist=c("GSM4390689","GSM4390690","GSM4390691")

for(i in 1:length(samplelist)){
samples=samplelist[i]
cell.file=paste0("/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/SRA/arab/",samples)

#when sample data is mouse_int
#cell.file=paste0("/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/SRA/mouse_int")

celltype=read.csv(paste0(cell.file,"/cell.csv"))

write.table(celltype, file = paste0(cell.file,"/","clusters_filterCB.txt"), 
          quote = FALSE, row.names = FALSE,col.names=FALSE,sep="\t")

annotation.files=paste0(cell.file,"/","clusters_filterCB.txt")



##add_cell.ident-------------------------------------------
#whether add cell.ident(sample.id) on cell barcode,if single sample,recommend FALSE; if merge peak from multiple samples,recommend TRUE
add_cell.ident=FALSE

##chr.modify------------------------------------------

###outfile.path----------------------------------------------
#outpath=paste0("/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/PAC/mouse_esc/scapa/",samples)
outpath=paste0("/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/PAC/arab/scapa/",samples)
temp.file=paste0(outpath,"/temp")
setwd(temp.file)

##mergebam---------------------------------------------------
merged.bam=paste0("dedup.",samples,".filterCB.bam")

##merge peaks------------------------------------------------------
#merge peaks less than 100 nt apart


scAPA::merge_peaks(bedtools.path = "", peaks.file = "Peakfile", path = "./")

##when sample data is arab can not run this
#***replace chr of merge.peakfile.bed,because scAPA mm10 (1-19,X,Y.MT)and h19 (1-22,X,Y,M)bed chr column no chr(1-19,X,Y.MT)*****######
#merge.peaks=read.delim(paste0(temp.file,'/merge.peakfile.bed') , header = F)
#merge.peaks$V1=gsub("chr","",merge.peaks$V1)
#scAPA::write.bed(.x = merge.peaks, f = "nochr.merge.peaks.bed")

#intersect Peakfile with a 3' UTR bed file to create a bed file of the 3’ UTR peaks. The peaks are annotated according to their 3’ UTR and their position within it
                                    
 peaks.bed <- intersect_peaks2(org = "Arab", utr_bed=utr,
                                    bed.name = "merge.peakfile.bed",
                                    path = "", bedtools.path = "", 
                                    introns = F)                                   

##add chr to peaks.bed of chrom column, because merged.wig of chrom column contain chr(when sample is arab ,not run)

#head(peaks.bed$V1)
#peaks.bed$V1=paste0("chr",peaks.bed$V1)

scAPA::write.bed(.x = peaks.bed, f = "peaks.bed")

#Separating Peaks with bimodal UMI counts distribution-------------------

wig.plus.command <- paste0("bedtools genomecov -ibam ",merged.bam,
                           " -bg -strand + | awk 'BEGIN ",
                           "{OFS = \"\t\"}{print $1",
                           ", $2, $3, $4, \".\", \"+\"}' > merged.wig")
system(command = wig.plus.command, wait = T)

wig.minus.command <- paste0( "bedtools genomecov -ibam ",merged.bam ,
                            " -bg -strand - | awk 'BEGIN {OFS = \"\t\"}{print $1",
                            ", $2, $3, $4, \".\", \"-\"}' >> merged.wig")
system(command = wig.minus.command, wait = T)

#intersect them with the peaks' BED file---------------------------------------

intersect.wig.command <- paste0( "bedtools intersect -s -wb ",
                                "-b peaks.bed -a merged.wig > intersected.wig")
system(intersect.wig.command, wait = T)


##Separating Peaks--------------------------------------------------------------

peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)
bed <- plyr::rbind.fill(parallel::mclapply(1:length(peaks.wig),
                                           FUN = creat_mclus,
                                           mc.cores = 1,
                                           mc.preschedule = T))

utr.saf <- bed[, c(4, 1, 2, 3, 6)]
rm("bed")
colnames(utr.saf) <- c("GeneID", "Chr", "Start", "End", "Strand")


# 2.Quantifying the usage of each peak -------------
                            
##featureCounts for 3UTRs for each cell------------------------------------
 if (!dir.exists("CellBams"))
        dir.create("CellBams")
   
    full.list.cells <- data.frame() 
                                      
  for (j in 1:length(annotation.files)) {
          sample=samples[j]
 
        list.cells <- read.delim(file = annotation.files[j], header = F)
        list.cells[, 1] <- as.character(list.cells[, 1])
        list.cells <- as.data.frame(list.cells)
        list.cells$Cell <- paste0(sample, "_", list.cells$V1)
        full.list.cells <- rbind.data.frame(full.list.cells,
                                            list.cells[, c(3, 2)])
        list.cells <- split(x = list.cells, f = list.cells$V1, drop = T)
        split_bams <- function(x) {
            cell <- x[, 1]
            bam <- paste0(x[, 3], ".bam")
            FilterBAMbyTag.command <- paste0(
                                             "/mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/softwares/Drop-seq_tools-2.2.0/FilterBamByTag ",
                                             "TAG=CB I=",
                                             merged.bam, " O=",
                                             "CellBams/", bam,
                                             " TAG_VALUE=", cell)
            system(command = FilterBAMbyTag.command, wait = T)
        }

  }
  
  parallel::mclapply(X = list.cells, FUN = split_bams, mc.cores = 12,
                           mc.preschedule = T) 
                           
 if(add_cell.ident==TRUE){                          
##obtain each cell name 
bam.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                    full.names = F)
cellnames <- gsub(x = bam.files, pattern = ".bam", replacement = "")
}

 if(add_cell.ident==FALSE){                          
##obtain each cell name 
bam.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                    full.names = F)
cellnames <- gsub(x = bam.files, pattern = ".bam", replacement = "")
cellnames <- gsub(x=cellnames,pattern =paste0(sample, "_") , replacement = "")

full.list.cells$Cell <- gsub(x=full.list.cells$Cell,pattern =paste0(sample, "_") , replacement = "")
}
##obtain bam file path
bam.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                    full.names = T)
                                                                        
# Count reads ----------                                    
#largestOverlap = True is specified so that reads spanning two peaks are counted according to their largest overlap
counts <- Rsubread::featureCounts(files = bam.files, isGTFAnnotationFile = F,
                                  strandSpecific = 1, annot.ext = utr.saf,
                                  largestOverlap = T,nthreads = 12)

co <- cbind.data.frame(rownames(counts$counts), counts$counts)
colnames(co) <- c("Peak_ID", cellnames)

#produce an object storing the 200 nt downstream sequences for each peak
#chr.modify:Whether add "chr" in the chromosome names of utr.saf
aseq <- read_down.seq2(saf = utr.saf, 
                             char.length.path = char.length.path,
                             fasta.path = fasta.path, 
                             bedtools.path="",chr.modify = F)
aseq <- aseq[,c(4, 6)]


##creat a scAPAList object for downstream analysis
row.Data <- counts$annotation[,c(2,3,4,1,6,5)]

a <- scAPA::set_scAPAList(.cells.counts = co, .row.Data = row.Data, 
                          .down.seq = aseq,
                          .cluster.anot = full.list.cells)
                          
setwd(outpath)
if (!dir.exists("outs"))
        dir.create("outs")

write.table(utr.saf, file = "outs/out_peaks.txt", 
          quote = FALSE, row.names = FALSE,col.names=TRUE,sep="\t")
          
write.table(co, file = "outs/peaks_counts.txt", 
          quote = FALSE, row.names = FALSE,col.names=TRUE,sep="\t")         
          
saveRDS(object = a, file = "outs/cell_Peaks.RDS")
  } 
  
                                   
end.time=Sys.time()
print(end.time-start.time)                                                                                              

