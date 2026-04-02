#!/bin/bash

###Preprocessing utr and intron region for inferring alternative polyadenylation events

###mouse:genecodev25
bedtools sort -i /mnt/user/BXY/ref_file/GRcm38fa+gtf/gencode.vM25.annotation.gtf  | bgzip > /mnt/user/BXY/ref_file/GRcm38fa+gtf/gencode.vM25.annotation.gtf.gz

cd /mnt/user/BXY/ref_file/GRcm38fa+gtf/SCAPE_bed
python /mnt/user/BXY/scAPA/scAPA_tools/SCAPE-main/main.py prepare --gtf /mnt/user/BXY/ref_file/GRcm38fa+gtf/gencode.vM25.annotation.gtf.gz --prefix GRCm38_v25


##human:genecodev34
bedtools sort -i /mnt/user/BXY/ref_file/GRch38fa+gtf/genecode/v34/gencode.v34.annotation.gtf  | bgzip > /mnt/user/BXY/ref_file/GRch38fa+gtf/genecode/v34/gencode.v34.annotation.gtf.gz

cd /mnt/user/BXY/ref_file/GRch38fa+gtf/genecode/v34
mkdir SCAPE_bed
cd SCAPE_bed
python /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/scAPA_tools/SCAPE-main/main.py prepare --gtf /mnt/user/BXY/ref_file/GRch38fa+gtf/genecode/v34/gencode.v34.annotation.gtf.gz --prefix GRCh38_v34


###Arabidopsis:tair10
bedtools sort -i /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/Arabidopsis_thaliana.TAIR10.49.gtf  | bgzip > /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/Arabidopsis_thaliana.TAIR10.49.gtf.gz

cd /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10
mkdir SCAPE_bed
cd SCAPE_bed
python /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/scAPA_tools/SCAPE-main/main.py prepare --gtf /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/Arabidopsis_thaliana.TAIR10.49.gtf.gz --prefix tair10


