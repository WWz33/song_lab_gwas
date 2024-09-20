#!/bin/bash 

mpheno="$1"
pheno_file="$2"
bfile="$3"
thread="$4"

## fam_file
## ID FID phe1 phe2 phe3
# 方便进行多表型分析，将其重命名为列名


# 读取文件的第一行（列名）
header=$(head -n 1 "$pheno_file")

# 使用空格分隔符将列名存储到数组中
IFS=' ' read -r -a columns <<< "$header"

# 根据 mpheno 参数选择相应的列名
if [ "$mpheno" -eq 1 ]; then
  out_var=${columns[2]}
else
  out_var=${columns[$((mpheno + 1))]}
fi

# 打印选择的列名
echo "Selected column for --out: $out_var"


# 生成新的 pheno 文件，去掉第一行并选择相应的列
new_pheno_file="${out_var}_pheno.txt"
awk -v mpheno="$mpheno" 'NR == 1 {next} {print $1, $2, $(mpheno + 2)}' "$pheno_file" > "$new_pheno_file"

head $new_pheno_file

gcta64 --mlma --bfile $bfile --pheno "$new_pheno_file" --mpheno 1  --maf 0.05 --out "$out_var" --thread-num $thread