# 定义一个函数来执行 Fisher 精确检验
perform_fisher_test <- function(input_file, output_file) {
  # 读取输入文件
  mydata1 <- read.table(input_file, header = FALSE, sep = "\t")
  
  # 获取行数
  size1 <- dim(mydata1)[1]
  
  # 打开输出文件以进行写入
  con <- file(output_file, "w")
  
  # 循环遍历每一行
  for (i in 1:size1) {
    # 提取当前行的值
    mC_HT   <- mydata1[i, 4]
    mC_NT   <- mydata1[i, 5]
    totalC_HT <- mydata1[i, 6]
    totalC_NT <- mydata1[i, 7]
    
    # 计算未甲基化值
    umC_HT <- totalC_HT - mC_HT
    umC_NT <- totalC_NT - mC_NT
    
    # 创建 2x2 矩阵用于 Fisher 精确检验
    data <- matrix(c(mC_HT, mC_NT, umC_HT, umC_NT), byrow = TRUE, ncol = 2)
    
    # 执行 Fisher 精确检验
    testResult <- chisq.test(data)
    
    # 提取 p 值
    Pvalue <- testResult$p.value
    
    # 创建输出数据
    outdata <- paste(mydata1[i, 1], mydata1[i, 2], mydata1[i, 3], Pvalue, sep = "\t")
    
    # 将输出数据写入文件
    writeLines(outdata, con)
  }
  
  # 关闭输出文件
  close(con)
}

# 从命令行参数获取输入和输出路径
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) != 2) {
  stop("Usage: Rscript script.R <input_file> <output_file>")
}

# 获取输入和输出文件路径
input_path <- args[1]
output_path <- args[2]

# 调用函数
perform_fisher_test(input_path, output_path)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript script.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

mydata1 <- read.table(input_file, header = FALSE, sep = "\t")

p <- mydata1$V4  

myresult <- p.adjust(p, method = "fdr")

mydata1$V5 <- myresult  

write.table(mydata1[myresult < 0.05, ], file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

