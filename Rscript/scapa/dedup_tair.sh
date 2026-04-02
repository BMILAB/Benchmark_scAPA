######Record runtime and memory usage
#/usr/bin/time -f "\n===== scAPA性能数据(Dedup Bam) =====\n峰值内存: %M KB\n平均内存: %t KB\n平均总内存(data+stack+text): %K KB\nCPU使用率: %P\n运行时长: %e秒\n===================" \
#    -a -o /mnt/T20T-6d12/BXY/scAPA/time.test/mouse.sperm/scapa/log.txt bash /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/code/scAPA/dedup.sh


#!/bin/bash
##output path

for sample in pas1_gn1000_bn1400_rep1 pas1_gn1000_bn1400_rep2 pas1_gn1000_bn1400_rep3; 
do echo $sample;
cd /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/tair_GSM3490690/scapa
mkdir ${sample}
mkdir ${sample}/temp
cd ${sample}/temp
##filterbam and PCR dedup
/mnt/user/BXY/software/Drop-seq_tools-2.2.0/FilterBam TAG_RETAIN=UB I=/mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/tair_GSM3490690/${sample}.bam O=UB.${sample}.bam;

samtools index UB.${sample}.bam;

umi_tools dedup -I UB.${sample}.bam -S dedup.${sample}.bam --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB;

##findpeaks
makeTagDirectory Tagdirectory dedup.*.bam; 
findPeaks Tagdirectory -size 50 -fragLength 115 -minDist 1 -strand separate -o Peakfile;
done

 

