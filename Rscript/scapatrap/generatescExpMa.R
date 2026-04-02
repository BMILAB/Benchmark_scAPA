
####************Aim*********************---------------------------------------------------------------------------------------
###This step combines peaks’ meta data (chr/strand/start/end) and peaks’ count data (peak-cell count matrix) to a scAPAtrapData. scAPAtrapData is stored in scAPAtrapData.rda, which is a list containing two matrices– peaks.meta and peaks.count.

##Here, we store the peaks' metadata (from the peak calling step) and peak count information in scAPAtrapData to generate APA data that has not undergone the pA anchoring step (i.e., extracting reads with A/T stretches to help locate precise pA sites).

###-------------------------------------------------------------------------------------------------------------------------



library(scAPAtrap)


## full path of tools-----
tools=list(samtools='/home/test/miniconda3/envs/R4.3/bin/samtools',
           umitools='/home/test/miniconda3/envs/R4.3/bin/umi_tools',
           featureCounts="/home/test/miniconda3/envs/R4.3/bin/featureCounts",
           star='/home/test/miniconda3/envs/R4.3/bin/STAR')
           

## input path--------------------------------------------------------------
###mouse_sperm
#dir0='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/mouse_sperm'

###human_covid_pbmc
#dir0='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/human_GSM4712885'

###tair
dir0='/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/tair_GSM3490690'

## output dir------------------------------------------------------------

###mouse_sperm
#out_path="/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/mouse_sperm/scapatrap"
#GSMlists=c("pas1_gn5000_bn2000_rep1","pas1_gn5000_bn2000_rep2","pas1_gn5000_bn2000_rep3")

###human_covid_pbmc
#out_path="/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/human_GSM4712885/scapatrap"
#GSMlists=c("pas1_gn5000_bn3000_rep1","pas1_gn5000_bn3000_rep2","pas1_gn5000_bn3000_rep3")

###tair_GSM3490690
out_path="/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/tair_GSM3490690/scapatrap"
GSMlists=c("pas1_gn1000_bn1400_rep1","pas1_gn1000_bn1400_rep2","pas1_gn1000_bn1400_rep3")   #



#read2 length----------------------------------------------------------------------------------------------------
#read_lens=c(98,98,98) ###mouse_sperm

#read_lens=c(91,91,91)  ###human_covid_pbmc

read_lens=c(115,115,115)  ###tair_GSM3490690


for(i in 1:length(GSMlists)){
GSM=GSMlists[i]

#inputBam=paste0(dir0,"/",GSM, "/filterCB.bam")
inputBam=paste0(dir0,"/",GSM, ".bam")

## log file to LOG all information (time, command, output file names..)
logf=gsub('.bam', '.APA.notails.onestep.log', inputBam, fixed=TRUE)


#output path will be under inputBam's dir if only dirname is provided(only provided name,not created); otherwise use the full path)
output_dir <- paste0(out_path,"/",GSM)
#output_dir <- paste0(out_path)

# print(output_dir)
if (!dir.exists(output_dir)){
   dir.create(output_dir)
} else {
    print("Dir already exists!")
} 

outputDir=paste0(output_dir,"/APA.notails")
#outputDir=paste0(output_dir)

if (!dir.exists(outputDir)){
   dir.create(outputDir)
} else {
    print("outputDir already exists!")
} 

##set parameters---------------------------------------------------------------------------------------------------------------
trap.params=setTrapParams()



#trap.params$barcode <- read.delim2(paste0(dir0, "/",GSM,'/barcode.txt'), header = F)$V1
trap.params$barcode <- read.delim2(paste0(dir0, "/",GSM,'.barcode_list.txt'), header = F)$V1

#delete "-1" of barcode
trap.params$barcode <-gsub("-[0-9]","",trap.params$barcode)

## set tail search way("genome","peaks","no")
trap.params$tails.search='no'

#chr.name
#trap.params$chrs=paste0("chr",c(as.character(1:19),'X','Y')) #mouse
#trap.params$chrs=paste0("chr",c(as.character(1:22),'X','Y'))  #human
trap.params$chrs=c(as.character(1:5))                        #tair


#bam contain CB,UB?
trap.params$TenX=TRUE

trap.params$maxwidth=1000

#read2 length
trap.params$readlength=read_lens[i]

trap.params$cov.cutoff=2

trap.params$min.cells=2

trap.params$min.count=2
trap.params$thread=12


## generatescExpMa: generate final PA data including peaks and counts information------------------

countsfile=paste0(output_dir,"/APA.tails/counts.tsv.gz")
peaksfile=paste0(output_dir,"/APA.tails/peaks.saf")
outputfile=paste0(outputDir, '/scAPAtrapData.rda')


outputfile <- generatescExpMa(countsfile=countsfile,
                              peaksfile=peaksfile,
                              barcode=trap.params$barcode,
                              #tails=tailsfile, 
                              #d=trap.params$d,
                              min.cells=trap.params$min.cells, 
                              min.count=trap.params$min.count,
                              ofile=outputfile, 
                              logf=logf)
## Final scAPAtrapData.rda
scAPAtrap:::.checkAndPrintFiles(varnames='outputfile', 
                                filenames=outputfile, logf=logf)
                                
}


