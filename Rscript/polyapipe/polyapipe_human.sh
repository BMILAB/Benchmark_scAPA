#!/bin/bash

# Set the log file path
LOG_FILE="/mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/human_GSM4712885/polyapipe/performance_log.txt"

# Creat the log file 
if [ ! -f "$LOG_FILE" ]; then
  echo "===== polyApipe 性能监控日志 =====" > "$LOG_FILE"
  echo "创建时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  echo "=============================" >> "$LOG_FILE"
fi


for sample in pas1_gn5000_bn3000_rep1 pas1_gn5000_bn3000_rep2 pas1_gn5000_bn3000_rep3;
do echo $sample;
mkdir /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/human_GSM4712885/polyapipe/${sample};
cd /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/human_GSM4712885/polyapipe/${sample}

echo -e "\n\n===== 新一轮分析开始 =====" >> "$LOG_FILE"
  echo "样本: $sample" >> "$LOG_FILE"
  echo "启动时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  
echo -e "\n[polyApipe.py 执行]" >> "$LOG_FILE"
  start_time=$(date +%s.%N)
 
  /usr/bin/time -f "峰值内存: %M KB\n平均内存: %K KB\nCPU使用率: %P\n运行时长: %e秒" \
    -a -o "$LOG_FILE" \
/mnt/user/BXY/scAPA/scAPA_tools/polyApipe-master/polyApipe.py -i /mnt/T20T-d2bf/BXY/scAPA/sim.data/SRA/human_GSM4712885/${sample}.bam -o /mnt/T20T-d2bf/BXY/scAPA/sim.data/PAC/human_GSM4712885/polyapipe/${sample}/${sample} --t 12 --cell_barcode_tag CB --umi_tag UB


end_time=$(date +%s.%N)
  elapsed=$(echo "scale=2; ($end_time - $start_time) / 60" | bc) 
  echo "实际运行时长: $elapsed 分钟" >> "$LOG_FILE"

  total_time=$(echo "$elapsed" | bc)
  max_peak=$(grep "峰值内存" "$LOG_FILE" | tail -1 | head -1 | awk '{print $2}')
  
  
  echo -e "\n===== 本轮汇总 =====" >> "$LOG_FILE"
  echo "总运行时长: $total_time 分钟" >> "$LOG_FILE"
  echo "最大峰值内存: $max_peak KB ($(echo "$max_peak/1024" | bc) MB)" >> "$LOG_FILE"
  echo "完成时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  echo "=============================" >> "$LOG_FILE"


done








