import csv

# 输入和输出文件名
input_file = 'my_gene_tpm_matrix.csv'
output_file = 'my_gene_tpm_matrix.txt'

# 打开输入文件和输出文件
with open(input_file, 'r') as csvfile, open(output_file, 'w', newline='') as txtfile:
    # 创建 CSV 读取器，指定逗号为分隔符
    reader = csv.reader(csvfile, delimiter=',')
    
    # 创建制表符分隔的写入器
    writer = csv.writer(txtfile, delimiter='\t')
    
    # 逐行读取并写入
    for row in reader:
        writer.writerow(row)

print(f"Converted {input_file} to {output_file} with tab separation.")
