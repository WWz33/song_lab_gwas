#!/bin/bash 
## for i in {1..6};do bash reml.sh $i alk_death_radio.phe  snp_removeIndividual 24 ;done
mpheno="$1"
pheno_file="$2"
grm="$3"
thread="$4"

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

gcta64 --grm $grm --pheno "$new_pheno_file"  --reml  --mpheno 1 --thread-num $4 --out "$out_var"

rm $new_pheno_file