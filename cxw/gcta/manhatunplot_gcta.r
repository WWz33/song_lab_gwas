#!/usr/bin/env Rscript

library(optparse)
library(qqman)
library(data.table)
# 定义命令行参数
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default = NULL,
              help = "Path to the directory containing GWAS files", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "Output directory", metavar = "character")
)
## e.g Rscript manhatunplot_gcta.r --path . --output_dir ./plot2
# 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 检查是否提供了必要的参数
if (is.null(opt$path) || is.null(opt$output_dir)) {
  print_help(opt_parser)
  stop("Both --path and --output_dir must be supplied", call. = FALSE)
}

# 获取输入路径和输出目录
input_path <- opt$path
output_dir <- opt$output_dir

gwas_list <- list.files(path = input_path, pattern = ".mlma$", full.names = TRUE)

output_dir <- output_dir

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 循环读取每个文件并生成 QQ 图和曼哈顿图
for (file in gwas_list) {
  # 读取数据
  gwas_data <- fread(file)
  
  # 提取文件名（不包括路径和扩展名）
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # 生成 QQ 图
  qq_file <- file.path(output_dir, paste0(file_name, "_qqplot.png"))
  png(qq_file)
  qq(gwas_data$p,main = file_name )
  dev.off()
  
  # 生成曼哈顿图
  gwas_data <- gwas_data[!(is.na(gwas_data$p)),]
  manhattan_file <- file.path(output_dir, paste0(file_name, "_manhattan.pdf"))
  pdf(manhattan_file,height = 6,width = 12)
  manhattan(gwas_data,chr = "Chr",bp = "bp",p = "p",snp = "SNP",logp = T, main=file_name)
  dev.off()
}
