#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip

def gene_set_to_bed(args):
    print('making gene set bed file')
    GeneSet = pd.read_csv(args.gene_set_file, header = None, names = ['GENE'])
    all_genes = pd.read_csv(args.gene_coord_file, sep='\s+')
    df = pd.merge(GeneSet, all_genes, on = 'GENE', how = 'inner')
    #shift the -1 to here
    df['START'] = np.maximum(1, df['START'] - args.windowsize) - 1
    df['END'] = df['END'] + args.windowsize
    df=df.drop(columns=['GENE'])
    pr_data=df.sort_values(by=['CHR', 'START', 'END'])
    merged_stack=[]
    for chro in pr_data['CHR'].unique():
        df=pr_data.loc[pr_data['CHR']==chro]
        stack=[]
        intervals=df[['START','END']].values.tolist()
        stack.append(intervals[0])
        for i in intervals[1:]:
            if stack[-1][0] <= i[0] <= stack[-1][-1]:
                stack[-1][-1] = max(stack[-1][-1], i[-1])
            else:
                stack.append(i)
    
        for e in stack:
            e.insert(0,chro)
        merged_stack = merged_stack + stack
    bed_for_annot = pd.DataFrame(merged_stack,columns=['CHR','START','END'])
    return bed_for_annot

def make_annot_files(args, bed_for_annot):
    print('making annot file')
    df_bim = pd.read_csv(args.bimfile, sep='\s+', usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = pd.DataFrame([['chr'+str(x1), x2 - 1, x2] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])])
    merged_stack=[]
    for chro in bed_for_annot['CHR'].unique():
        df_a=bed_for_annot.loc[bed_for_annot['CHR']==chro]
        intervals_a=df_a[['START','END']].values.tolist()
        df_b = iter_bim[iter_bim[0]==chro]
        intervals_b=df_b[[1,2]].values.tolist()
        i = j = 0
         
        n = len(intervals_a)
        m = len(intervals_b)
        stack=[]
        while i < n and j < m:
            # Left bound for intersecting segment
            l = max(intervals_a[i][0], intervals_b[j][0])
            # Right bound for intersecting segment
            r = min(intervals_a[i][1], intervals_b[j][1])
            # If segment is valid print it
            if l <= r: 
                stack.append([l,r])
                     
            # If i-th interval's right bound is 
            # smaller increment i else increment j
            if intervals_a[i][1] < intervals_b[j][1]:
                i += 1
            else:
                j += 1

        for e in stack:
            e.insert(0,chro)
        merged_stack = merged_stack + stack

    annot_bed = pd.DataFrame(merged_stack,columns=['CHR','START','END'])
    df_int = pd.DataFrame({'CHR': annot_bed.CHR.str.lstrip('chr').astype(np.int64),'BP': annot_bed.END, 'ANNOT':1})
    df_annot = pd.merge(df_bim, df_int, how='left', on=['CHR','BP'])
    df_annot.fillna(0, inplace=True)
    df_annot = df_annot[['ANNOT']].astype(int)

    if args.annot_file.endswith('.gz'):
        with gzip.open(args.annot_file, 'wb') as f:
            df_annot.to_csv(f, sep = "\t", index = False)
    else:
        df_annot.to_csv(args.annot_file, sep="\t", index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-set-file', type=str, help='a file of gene names, one line per gene.')
    parser.add_argument('--gene-coord-file', type=str, default='ENSG_coord.txt', help='a file with columns GENE, CHR, START, and END, where START and END are base pair coordinates of TSS and TES. This file can contain more genes than are in the gene set. We provide ENSG_coord.txt as a default.')
    parser.add_argument('--windowsize', type=int, help='how many base pairs to add around the transcribed region to make the annotation?')
    parser.add_argument('--bed-file', type=str, help='the UCSC bed file with the regions that make up your annotation')
    parser.add_argument('--nomerge', action='store_true', default=False, help='don\'t merge the bed file; make an annot file wi    th values proportional to the number of intervals in the bedfile overlapping the SNP.')
    parser.add_argument('--bimfile', type=str, help='plink bim file for the dataset you will use to compute LD scores.')
    parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')

    args = parser.parse_args()

    if args.gene_set_file is not None:
        bed_for_annot = gene_set_to_bed(args)
    else:
        #### fix else clause below
        bed_for_annot = BedTool(args.bed_file).sort()
        if not args.nomerge:
            bed_for_annot = bed_for_annot.merge()

    make_annot_files(args, bed_for_annot)
