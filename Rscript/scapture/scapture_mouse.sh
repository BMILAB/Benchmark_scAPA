######Record the memory used and the time
#/usr/bin/time -f "\n=====SCAPTURE 性能数据 =====\n峰值内存: %M KB\n平均内存: %t KB\n平均总内存(data+stack+text): %K KB\nCPU使用率: %P\n运行时长: %e秒 \n===================" \
#    -a -o /mnt/user/BXY/scAPA/time.test/mouse.sperm/scapture/log.txt bash /mnt/user/BXY/scAPA/code/scapture/scapture.sh

#!/bin/bash

export PATH=/mnt/user/BXY/scapture/SCAPTURE-main:$PATH

### Run scapture annotation module

#INPUT file
#(1)genome fa and fai file
#(2)gene annotation file(GTF format) and chromosom size file
#output dir
cd /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/mouse_sperm/scapture


#	scapture -m annotation -o /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/SCAPTURE_anno -g /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/GRCm38.p6.genome.fa --gtf /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/gencode.vM25.annotation.gtf --cs /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/main.chrNameLength.txt --extend 2000 &> /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/SCAPTURE_anno/annotation.log
	
for sample in pas1_gn5000_bn2000_rep1 pas1_gn5000_bn2000_rep2 pas1_gn5000_bn2000_rep3;
do echo $sample;
mkdir /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/mouse_sperm/scapture/${sample};
cd /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/mouse_sperm/scapture/${sample}

### Run scapture PAScall module	

scapture -m PAScall -a /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/SCAPTURE_anno/SCAPTURE_annotation -g /mnt/user/BXY/ref_file/GRcm38fa+gtf/genecode/GRCm38.p6.genome.fa -b /mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/mouse_sperm/${sample}.bam  -l 98 -o ${sample} -p 12 --species mouse --polyaDB /mnt/user/BXY/ref_file/refPAC/mm10.SupTab_KnownPASs_fourDBs.txt &> ${sample}.PAScall.log;


### Run scapture PASquant module
# 1. Select PASs with positive prediction and konwn sites overlapped (positive prediction or known pA,recomanded)
	perl -alne '$,="\t";print @F[0..11] if $F[12] > 0 | $F[13] eq "positive";' ${sample}.3primeExtended.peaks.evaluated.bed ${sample}.exonic.peaks.evaluated.bed ${sample}.intronic.peaks.evaluated.bed > ${sample}_ALL.PASquant.bed;	
	

scapture -m PASquant -b /mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/mouse_sperm/${sample}.bam --celllist /mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/mouse_sperm/${sample}.barcode_list.txt --pas ${sample}_ALL.PASquant.bed -o ${sample}.PASquant &> ${sample}.PASquant.log;


done


