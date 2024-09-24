#!/usr/bin/python3

import argparse
import sys,re,os
import pandas as pd

def parse_arguments():
	parser = argparse.ArgumentParser(description='Extract candidate gene from GWAS output files')
	
	parser.add_argument('--input', help='GWAS output file', required=True, action="store")
	parser.add_argument('-b', help='Gene position (.bed)', required=True, action="store")
	parser.add_argument('-f', help='Protein sequence (.fasta)', required=True, action="store")
	parser.add_argument('-p', type=float, default=1e-5, help='Significant P-value', required=True, action="store")
	parser.add_argument('-d', type=int, default=100000, help='The distance in bp from a variant site to search for genes (default: 100000)', required=False, action="store")
	parser.add_argument('-n', type=int, default=3, help='The number of signals in QTL (default: 3)', required=False, action="store")
	parser.add_argument('-o', help='The prefix of output files.', required=True, action="store")
	
	return parser.parse_args()


def save_gwas_results(input_file,pvalue,distance):
	signal_info = {}
	signal_bed = {}
	with open(input_file) as f:
		for line in f.read().splitlines()[1:]:
			gwas_results = line.split('\t')
			var_id = gwas_results[1]
			if float(gwas_results[12]) <= pvalue:
				bed_region_start = int(gwas_results[2])-distance+1
				bed_region_end = int(gwas_results[2])+distance-1
				bed_region = [gwas_results[0], gwas_results[2], gwas_results[2],gwas_results[7],gwas_results[12]]
				signal_info[var_id] = gwas_results
				signal_bed[var_id] = '\t'.join(bed_region)
		
	return signal_info,signal_bed



def getsequence(fasta,bedfile):
        df = pd.read_csv(bedfile, delimiter='\t', header = None)
        geneid = df.iloc[:,8].tolist()
        seq={}
        gid=''
        for line in open(fasta):
                if line.startswith('>'):
                        line = line.replace('\n', '')
                        if '.cds' in line:
                                gid = line.replace('>','').split('.c')[0]
                        elif '.cds' not in line:
                                gid = line.replace('>','').split(' ')[0]
                        seq[gid] = ''
                else:
                        seq[gid] += line.replace('\n','')
        return seq,geneid


def output_sig_signal_bed(bed_dict, output):
	signal_bed_output = output+'.sig_gwas_signal.bed'
	fw_signal_bed_output = open(signal_bed_output, "w+")
	for k in bed_dict.keys():
		print(bed_dict[k],file = fw_signal_bed_output)
	return signal_bed_output


def overalaps_variant_gene(signal_bed_output,num,gene_file,pep_file,output):
	merged_bed_output = output+'.sig_gwas_QTL.bed'
	overlap_bed_output = output+'.sig_gwas_QTL_candidate.bed'
	candidate_fasta_output = open(output + '.sig_gwas_QTL_candidate.fasta', 'w+')
	command = "{print$0}"
	os.system("bedtools merge -i {0} -d 100000 -c 2,5 -o collapse,min | awk -F ',' 'NF >= {1} {2}' > {3}".format(signal_bed_output,num,command,merged_bed_output))
	os.system("bedtools intersect -a {0} -b {1} -wa -wb > {2}".format(merged_bed_output,gene_file,overlap_bed_output))
	f,b = getsequence(pep_file,overlap_bed_output)
	for k in b:
		gene = k+".t1"
		candidate_seq = '>'+gene+'\n'+f[gene]+'\n'
		candidate_fasta_output.write(candidate_seq)
	
	
if __name__ == "__main__":
	args = parse_arguments()
	input_file = args.input
	gene_file = args.b
	pep_file = args.f
	pvalue = args.p
	distance = args.d
	num = args.n
	output = args.o

	signal_dict, bed_dict = save_gwas_results(input_file,pvalue,distance)
	overalaps_variant_gene(output_sig_signal_bed(bed_dict, output),num,gene_file,pep_file,output)

