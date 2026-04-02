#!/bin/bash

# Set the log file path
LOG_FILE="/mnt/T20T-6d12/BXY/scAPA/sim_data/PAC/tair_GSM3490690/scape/performance_log.txt"

# Creat the log file
if [ ! -f "$LOG_FILE" ]; then
  echo "===== SCAPE 性能监控日志 =====" > "$LOG_FILE"
  echo "创建时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  echo "=============================" >> "$LOG_FILE"
fi


###apamix-----
##Perform the actual alternative polyadenylation events inference. If one wants to change the default paramters such as insert size distribution, maximum number of pA sites etc, it can be done by modifying the file "apamix.py" under folder "apamix".

##output folder not need create in advance
cd /mnt/T20T-6d12/BXY/scAPA/sim_data/SRA/tair_GSM3490690
for sample in pas1_gn1000_bn1400_rep1 pas1_gn1000_bn1400_rep2 pas1_gn1000_bn1400_rep3 ;
do echo $sample;
 
  echo -e "\n\n===== 新一轮分析开始 =====" >> "$LOG_FILE"
  echo "样本: $sample" >> "$LOG_FILE"
  echo "启动时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  
  echo -e "\n[main.py 执行]" >> "$LOG_FILE"
  start_time_1=$(date +%s.%N)
  
  /usr/bin/time -f "峰值内存: %M KB\n平均内存: %K KB\nCPU使用率: %P\n运行时长: %e秒" \
    -a -o "$LOG_FILE" \
python /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/scAPA_tools/SCAPE-main/main.py apamix --bed /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/ref_file/TAIR10/SCAPE_bed/tair10_utr.bed --bam ${sample}.bam --out /mnt/T20T-6d12/BXY/scAPA/sim_data/PAC/tair_GSM3490690/scape/${sample} --cores 12 --cb /mnt/T20T-6d12/BXY/scAPA/sim_data/SRA/tair_GSM3490690/${sample}.barcode_list.txt

end_time_1=$(date +%s.%N)
  elapsed_1=$(echo "scale=2; ($end_time_1 - $start_time_1)/60" | bc)  
  echo "实际运行时长: $elapsed_1 分钟" >> "$LOG_FILE"

  echo -e "\n[group_pa.py 执行]" >> "$LOG_FILE"
  start_time_2=$(date +%s.%N)
  
  /usr/bin/time -f "峰值内存: %M KB\n平均内存: %K KB\nCPU使用率: %P\n运行时长: %e秒" \
    -a -o "$LOG_FILE" \
python /mnt/64cf3476-350c-46ad-bc48-574fa64a0334/test/BXY/scAPA/scAPA_tools/SCAPE-main/scripts/group_pa.py --files /mnt/T20T-6d12/BXY/scAPA/sim_data/PAC/tair_GSM3490690/scape/${sample}/pasite.csv.gz --labels ${sample} --outfile /mnt/T20T-6d12/BXY/scAPA/sim_data/PAC/tair_GSM3490690/scape/${sample}/collapse_pa.tsv.gz --group_threshold 100


 end_time_2=$(date +%s.%N)
 elapsed_2=$(echo "scale=2; ($end_time_2 - $start_time_2)/60" | bc)
  echo "实际运行时长: $elapsed_2 分钟" >> "$LOG_FILE"
  
##Summary
  total_time=$(echo "$elapsed_1 + $elapsed_2" | bc)
  peak_mem_1=$(grep "峰值内存" "$LOG_FILE" | tail -2 | head -1 | awk '{print $2}')
  peak_mem_2=$(grep "峰值内存" "$LOG_FILE" | tail -1 | awk '{print $2}')
  max_peak=$(echo "$peak_mem_1 $peak_mem_2" | awk '{if($1>$2) print $1; else print $2}')
  
  echo -e "\n===== 本轮汇总 =====" >> "$LOG_FILE"
  echo "总运行时长: $total_time 分钟" >> "$LOG_FILE"
  echo "最大峰值内存: $max_peak KB ($(echo "$max_peak/1024" | bc) MB)" >> "$LOG_FILE"
  echo "完成时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  echo "=============================" >> "$LOG_FILE"
done

echo -e "\n监控结果已追加到: $LOG_FILE"

##Generate a report
echo "===== 执行报告 ====="
echo "样本: $sample"
echo "总运行时长: $total_time 分钟"
echo "最大峰值内存: $max_peak KB ($(echo "$max_peak/1024" | bc) MB)"
echo "完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "==================="




