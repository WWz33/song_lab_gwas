#!/bin/bash

# 显示帮助信息的函数
show_help() {
    echo "Usage: ld.sh -p <pre> -f <format> -w <win> -s <step> -r <r2>"
    echo "Options:"
    echo "  -p <pre>       Prefix for input and output files"
    echo "  -f <format>    Format of the input file (vcf, gz, bfile, plink)"
    echo "  -w <win>       Window size for LD pruning"
    echo "  -s <step>      Step size for LD pruning"
    echo "  -r <r2>        r² threshold for LD pruning"
    echo "  -h             Show this help message"
}

# 初始化变量
pre=""
format=""
win=""
step=""
r2=""

# 解析命令行选项
while getopts ":p:f:w:s:r:h" option; do
    case $option in
        p) pre=$OPTARG ;;
        f) format=$OPTARG ;;
        w) win=$OPTARG ;;
        s) step=$OPTARG ;;
        r) r2=$OPTARG ;;
        h) show_help
           exit 0 ;;
        \?) echo "Invalid option: -$OPTARG" >&2
            show_help
            exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2
           show_help
           exit 1 ;;
    esac
done

# 检查是否提供了必要的参数
if [ -z "$pre" ] || [ -z "$format" ] || [ -z "$win" ] || [ -z "$step" ] || [ -z "$r2" ]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# 根据格式执行相应的命令
case $format in
    vcf)
        plink --vcf ${pre}.vcf --indep-pairwise $win $step $r2 --allow-extra-chr --allow-no-sex --make-bed --out $pre
        plink --vcf ${pre}.vcf --extract ${pre}.prune.in --recode --allow-extra-chr --allow-no-sex --make-bed --out ${pre}_ldpruned
        ;;
    gz)
        plink --vcf ${pre}.vcf.gz --indep-pairwise $win $step $r2 --allow-extra-chr --allow-no-sex --make-bed --out $pre
        plink --vcf ${pre}.vcf.gz --extract ${pre}.prune.in --recode --allow-extra-chr --allow-no-sex --make-bed --out ${pre}_ldpruned
        ;;
    bfile)
        plink --bfile ${pre} --indep-pairwise $win $step $r2 --allow-extra-chr --allow-no-sex --make-bed --out $pre
        plink --bfile ${pre} --extract ${pre}.prune.in --recode --allow-extra-chr --allow-no-sex --make-bed --out ${pre}_ldpruned
        ;;
    plink)
        plink --file ${pre} --indep-pairwise $win $step $r2 --allow-extra-chr --allow-no-sex --make-bed --out $pre
        plink --file ${pre} --extract ${pre}.prune.in --recode --allow-extra-chr --allow-no-sex --make-bed --out ${pre}_ldpruned    
        ;;
    *)
        echo "Unsupported file type: $format"
        show_help
        exit 1
        ;;
esac