# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 10:22:47 2024

@author: antsh
"""
import sys, argparse
print(sys.path)
sys.path.append('C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\git\\ldsc_py3\\ldsc_py3')


import ldscore.ldscore as ld
import ldscore.parse as ps
import numpy as np
import pandas as pd
from subprocess import call
import gzip

#from ldscore.ldscore import block_left_to_right

args = parser.parse_args('')
##python ldsc.py --bfile /nas/longleaf/home/abrantes/ldsc_py3R/01data/toy_ldsc/1000G.EUR.QC.22 --l2 --ld-wind-cm 1 
##--out /nas/longleaf/home/abrantes/ldsc_py3R/03output/ldsc/toy_data/1000G.EUR.QC.22

args.bfile="C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\toy_ldsc\\1000G_Phase3_plinkfiles\\1000G_EUR_Phase3_plink\\1000G.EUR.QC.22"
args.out="C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\03output\\1000G.EUR.QC.22"
args.l2=True
args.ld_wind_cm=1

if args.bfile:
    snp_file, snp_obj = args.bfile+'.bim', ps.PlinkBIMFile
    ind_file, ind_obj = args.bfile+'.fam', ps.PlinkFAMFile
    array_file, array_obj = args.bfile+'.bed', ld.PlinkBEDFile
    
array_snps = snp_obj(snp_file)
m = len(array_snps.IDList)
    
annot_matrix, annot_colnames, keep_snps = None, None, None,
n_annot = 1

# read fam
array_indivs = ind_obj(ind_file)
n = len(array_indivs.IDList)
    
keep_indivs = None
    
geno_array = array_obj(array_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=args.maf)
 
x = np.array((args.ld_wind_snps, args.ld_wind_kb, args.ld_wind_cm), dtype=bool)

max_dist = args.ld_wind_cm
coords = np.array(array_snps.df['CM'])[geno_array.kept_snps]

block_left = ld.getBlockLefts(coords, max_dist)

scale_suffix = ''

lN = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_matrix)
col_prefix = "L2"; file_suffix = "l2"

ldscore_colnames = [col_prefix+scale_suffix]

out_fname = args.out + '.' + file_suffix + '.ldscore'
new_colnames = geno_array.colnames + ldscore_colnames
df = pd.DataFrame.from_records(np.c_[geno_array.df, lN])
df.columns = new_colnames
print(out_fname)

if args.print_snps:
    if args.print_snps.endswith('gz'):
        print_snps = pd.read_csv(args.print_snps, header=None, compression='gzip')
    elif args.print_snps.endswith('bz2'):
        print_snps = pd.read_csv(args.print_snps, header=None, compression='bz2')
    else:
        print_snps = pd.read_csv(args.print_snps, header=None)
    if len(print_snps.columns) > 1:
        raise ValueError('--print-snps must refer to a file with a one column of SNP IDs.')
    #log.log('Reading list of {N} SNPs for which to print LD Scores from {F}'.format(\
     #               F=args.print_snps, N=len(print_snps)))

    print_snps.columns=['SNP']
    df = df.ix[df.SNP.isin(print_snps.SNP),:]
    if len(df) == 0:
        raise ValueError('After merging with --print-snps, no SNPs remain.')
    else:
        msg = 'After merging with --print-snps, LD Scores for {N} SNPs will be printed.'
        #log.log(msg.format(N=len(df)))

l2_suffix = '.gz'
#log.log("Writing LD Scores for {N} SNPs to {f}.gz".format(f=out_fname, N=len(df)))
df.drop(['CM','MAF'], axis=1).to_csv(out_fname, sep="\t", header=True, index=False,
    float_format='%.3f')
#call(['gzip', '-f', out_fname])

with open(out_fname, 'rb') as f_in, gzip.open(out_fname+'.gz', 'wb') as f_out:
    f_out.writelines(f_in)
    
M = [geno_array.m]
M_5_50 = [np.sum(geno_array.maf > 0.05)]

print(df.head())
print(M)
print(M_5_50)
# print .M
#py3
with open(args.out + '.' + file_suffix + '.M', 'wb') as fout_M:
    ## added list to py3 changed functionality of map
    fout_M.write('\t'.join(list(map(str, M))).encode())
#py2
#fout_M = open(args.out + '.'+ file_suffix +'.M','wb')
#print >>fout_M, '\t'.join(map(str,M))
#fout_M.close()

# print .M_5_50
#py3
with open(args.out + '.' + file_suffix + '.M_5_50', 'wb') as fout_M_5_50:
    ## added list to py3 changed functionality of map
    fout_M_5_50.write('\t'.join(list(map(str, M_5_50))).encode())
#py2
#fout_M_5_50 = open(args.out + '.'+ file_suffix +'.M_5_50','wb')
#print >>fout_M_5_50, '\t'.join(map(str,M_5_50))
#fout_M_5_50.close()