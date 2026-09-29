#!/bin/bash  
  
# 定义路径  
input_dir1="/public/home/phwan/RNA-seq/7.5N_9NdiffExp/DESeq2.141481.dir/Result"  
input_dir2="/public/home/phwan/RNA-seq/6.HT_NT1diffExp/DESeq2.96681.dir/Result"  
output_dir="/public/home/phwan/RNA-seq/8.noback_DEG"  

for i in 6-7 7-8 8-9 9-10  
do
# 定义文件名  
base_file="59N${i}filter_result.txt"  
compare_file1="5HN${i}filter_result.txt"  
compare_file2="9HN${i}filter_result.txt"  
  
# 创建输出目录（如果不存在）  
mkdir -p "$output_dir"  
  
# 提取filter_result.txt中的基因ID  
awk '{print $1}' "$input_dir1/$base_file" > "$output_dir/gene_ids.txt"  
  
# 从5HN6-7filter_result.txt中找出不在59N6-7filter_result.txt中的基因ID  
awk 'NR==FNR {ids[$1]; next} !($1 in ids)' "$output_dir/gene_ids.txt" "$input_dir2/$compare_file1" > "$output_dir/5HN${i}DEG.txt"  
  
# 从9HN6-7filter_result.txt中找出不在59N6-7filter_result.txt中的基因ID  
awk 'NR==FNR {ids[$1]; next} !($1 in ids)' "$output_dir/gene_ids.txt" "$input_dir2/$compare_file2" > "$output_dir/9HN${i}DEG.txt"  
  
# 删除临时文件  
rm "$output_dir/gene_ids.txt"  
done  
echo "处理完成，结果已保存到 $output_dir/5HN${i}DEG.txt 和 $output_dir/9HN${i}DEG.txt"

