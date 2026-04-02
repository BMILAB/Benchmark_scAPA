
#!/bin/bash

#output folder
cd /media/ubantu_root/16TB_21/BXY/scAPA/PAC/arab/scraps/GSM4390690/sim_data_res
#mkdir pas1_gn5000_bn2000_rep1 pas1_gn5000_bn2000_rep2 pas1_gn5000_bn2000_rep3

# Set the log file path
LOG_FILE="/media/ubantu_root/16TB_21/BXY/scAPA/PAC/arab/scraps/GSM4390690/sim_data_res/performance_log.txt"

# Creat the log file 
if [ ! -f "$LOG_FILE" ]; then
  echo "===== scraps 性能监控日志 =====" > "$LOG_FILE"
  echo "创建时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  echo "=============================" >> "$LOG_FILE"
fi


for sample in GSM4390690;
do echo $sample;
cd /media/ubantu_root/16TB_21/BXY/scAPA/PAC/arab/scraps/${sample}

echo -e "\n\n===== 新一轮分析开始 =====" >> "$LOG_FILE"
  echo "样本: $sample" >> "$LOG_FILE"
  echo "启动时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  
  start_time=$(date +%s.%N)

##code run

/usr/bin/time -f "峰值内存: %M KB\n平均内存: %K KB\nCPU使用率: %P\n运行时长: %e秒" \
    -o "$LOG_FILE" -a \
    snakemake --snakefile Snakefile --configfile config.yaml --resources total_impact=12 --keep-going --cores 12 >> "$LOG_FILE"


end_time=$(date +%s.%N)
  elapsed=$(echo "scale=2; ($end_time - $start_time) / 60" | bc) 
  echo "实际运行时长: $elapsed 分钟" >> "$LOG_FILE"

  ##Summary Report
  total_time=$(echo "$elapsed" | bc)
  max_peak=$(grep "峰值内存" "$LOG_FILE" | tail -1 | head -1 | awk '{print $2}')
  
  
  echo -e "\n===== 本轮汇总 =====" >> "$LOG_FILE"
  echo "总运行时长: $total_time 分钟" >> "$LOG_FILE"
  echo "最大峰值内存: $max_peak KB ($(echo "$max_peak/1024" | bc) MB)" >> "$LOG_FILE"
  echo "完成时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG_FILE"
  echo "=============================" >> "$LOG_FILE"

done
