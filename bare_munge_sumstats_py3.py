# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 09:20:43 2024
Add p-values and test munge_sumstats script
@author: antsh
"""

#from __future__ import division
import pandas as pd
import numpy as np
import os
import traceback
import gzip
import bz2
from scipy.stats import chi2
#from ldscore import sumstats
#from ldsc import MASTHEAD, Logger, sec_to_str
import time
np.seterr(invalid='ignore')
import gc
#import csv
import sys, argparse
import itertools as it
from functools import reduce

print(sys.path)
sys.path.append('C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\git\\ldsc_py3\\ldsc_py3')

os.chdir("C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\git\\ldsc_py3\\ldsc_py3")

iq_df = pd.read_parquet('../../../01data/iq.sumstats/raw_sumstats.parquet')

print(iq_df.head())
iq_df_22=iq_df[iq_df.CHR==22]
print(iq_df_22.head())

in_df = pd.read_parquet('../../../01data/insomnia.sumstats/raw_sumstats.parquet')

print(in_df.head())
in_df_22=in_df[in_df.CHR=='22']
print(in_df.CHR.value_counts())
print(in_df_22.head())
print(in_df_22.columns)
dir(in_df_22)

iq_df_22.to_csv('../../../01data/iq.sumstats/iq_df_22.txt', sep='\t', index=False)
in_df_22.to_csv('../../../01data/insomnia.sumstats/insomnia_df_22.txt', sep='\t', index=False)

iq_df_22=pd.read_csv('../../../01data/insomnia.sumstats/insomnia_df_22.txt')


# =============================================================================
# with open("../../../01data/iq.sumstats/iq_df_22.txt", "w", newline="\n") as f:
#     writer = csv.writer(f, delimiter="\t")
#     writer.writerow(in_df_22.columns)
#     writer.writerows(in_df_22)
# =============================================================================




del iq_df
del in_df
gc.collect()

parser = argparse.ArgumentParser()
parser.add_argument('--sumstats', default=None, type=str,
                    help="Input filename.")
parser.add_argument('--N', default=None, type=float,
                    help="Sample size If this option is not set, will try to infer the sample "
                    "size from the input file. If the input file contains a sample size "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--N-cas', default=None, type=float,
                    help="Number of cases. If this option is not set, will try to infer the number "
                    "of cases from the input file. If the input file contains a number of cases "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--N-con', default=None, type=float,
                    help="Number of controls. If this option is not set, will try to infer the number "
                    "of controls from the input file. If the input file contains a number of controls "
                    "column, and this flag is set, the argument to this flag has priority.")
parser.add_argument('--out', default=None, type=str,
                    help="Output filename prefix.")
parser.add_argument('--info-min', default=0.9, type=float,
                    help="Minimum INFO score.")
parser.add_argument('--maf-min', default=0.01, type=float,
                    help="Minimum MAF.")
parser.add_argument('--daner', default=False, action='store_true',
                    help="Use this flag to parse Stephan Ripke's daner* file format.")
parser.add_argument('--daner-n', default=False, action='store_true',
                    help="Use this flag to parse more recent daner* formatted files, which "
		    "include sample size column 'Nca' and 'Nco'.")
parser.add_argument('--no-alleles', default=False, action="store_true",
                    help="Don't require alleles. Useful if only unsigned summary statistics are available "
                    "and the goal is h2 / partitioned h2 estimation rather than rg estimation.")
parser.add_argument('--merge-alleles', default=None, type=str,
                    help="Same as --merge, except the file should have three columns: SNP, A1, A2, "
                    "and all alleles will be matched to the --merge-alleles file alleles.")
parser.add_argument('--n-min', default=None, type=float,
                    help='Minimum N (sample size). Default is (90th percentile N) / 2.')
parser.add_argument('--chunksize', default=5e6, type=int,
                    help='Chunksize.')

# optional args to specify column names
parser.add_argument('--snp', default=None, type=str,
                    help='Name of SNP column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-cas-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--N-con-col', default=None, type=str,
                    help='Name of N column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--a1', default=None, type=str,
                    help='Name of A1 column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--a2', default=None, type=str,
                    help='Name of A2 column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--p', default=None, type=str,
                    help='Name of p-value column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--frq', default=None, type=str,
                    help='Name of FRQ or MAF column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--signed-sumstats', default=None, type=str,
                    help='Name of signed sumstat column, comma null value (e.g., Z,0 or OR,1). NB: case insensitive.')
parser.add_argument('--info', default=None, type=str,
                    help='Name of INFO column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--info-list', default=None, type=str,
                    help='Comma-separated list of INFO columns. Will filter on the mean. NB: case insensitive.')
parser.add_argument('--nstudy', default=None, type=str,
                    help='Name of NSTUDY column (if not a name that ldsc understands). NB: case insensitive.')
parser.add_argument('--nstudy-min', default=None, type=float,
                    help='Minimum # of studies. Default is to remove everything below the max, unless there is an N column,'
                    ' in which case do nothing.')
parser.add_argument('--ignore', default=None, type=str,
                    help='Comma-separated list of column names to ignore.')
parser.add_argument('--a1-inc', default=False, action='store_true',
                    help='A1 is the increasing allele.')
parser.add_argument('--keep-maf', default=False, action='store_true',
                    help='Keep the MAF column (if one exists).')




null_values = {

    'LOG_ODDS': 0,
    'BETA': 0,
    'OR': 1,
    'Z': 0
}

default_cnames = {

    # RS NUMBER
    'SNP': 'SNP',
    'MARKERNAME': 'SNP',
    'SNPID': 'SNP',
    'RS': 'SNP',
    'RSID': 'SNP',
    'RS_NUMBER': 'SNP',
    'RS_NUMBERS': 'SNP',
    # NUMBER OF STUDIES
    'NSTUDY': 'NSTUDY',
    'N_STUDY': 'NSTUDY',
    'NSTUDIES': 'NSTUDY',
    'N_STUDIES': 'NSTUDY',
    # P-VALUE
    'P': 'P',
    'PVALUE': 'P',
    'P_VALUE':  'P',
    'PVAL': 'P',
    'P_VAL': 'P',
    'GC_PVALUE': 'P',
    # ALLELE 1
    'A1': 'A1',
    'ALLELE1': 'A1',
    'ALLELE_1': 'A1',
    'EFFECT_ALLELE': 'A1',
    'REFERENCE_ALLELE': 'A1',
    'INC_ALLELE': 'A1',
    'EA': 'A1',
    # ALLELE 2
    'A2': 'A2',
    'ALLELE2': 'A2',
    'ALLELE_2': 'A2',
    'OTHER_ALLELE': 'A2',
    'NON_EFFECT_ALLELE': 'A2',
    'DEC_ALLELE': 'A2',
    'NEA': 'A2',
    # N
    'N': 'N',
    'NCASE': 'N_CAS',
    'CASES_N': 'N_CAS',
    'N_CASE': 'N_CAS',
    'N_CASES': 'N_CAS',
    'N_CONTROLS': 'N_CON',
    'N_CAS': 'N_CAS',
    'N_CON': 'N_CON',
    'N_CASE': 'N_CAS',
    'NCONTROL': 'N_CON',
    'CONTROLS_N': 'N_CON',
    'N_CONTROL': 'N_CON',
    'WEIGHT': 'N',  # metal does this. possibly risky.
    # SIGNED STATISTICS
    'ZSCORE': 'Z',
    'Z-SCORE': 'Z',
    'GC_ZSCORE': 'Z',
    'Z': 'Z',
    'OR': 'OR',
    'B': 'BETA',
    'BETA': 'BETA',
    'LOG_ODDS': 'LOG_ODDS',
    'EFFECTS': 'BETA',
    'EFFECT': 'BETA',
    'SIGNED_SUMSTAT': 'SIGNED_SUMSTAT',
    # INFO
    'INFO': 'INFO',
    # MAF
    'EAF': 'FRQ',
    'FRQ': 'FRQ',
    'MAF': 'FRQ',
    'FRQ_U': 'FRQ',
    'F_U': 'FRQ',
}

describe_cname = {
    'SNP': 'Variant ID (e.g., rs number)',
    'P': 'p-Value',
    'A1': 'Allele 1, interpreted as ref allele for signed sumstat.',
    'A2': 'Allele 2, interpreted as non-ref allele for signed sumstat.',
    'N': 'Sample size',
    'N_CAS': 'Number of cases',
    'N_CON': 'Number of controls',
    'Z': 'Z-score (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'OR': 'Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)',
    'BETA': '[linear/logistic] regression coefficient (0 --> no effect; above 0 --> A1 is trait/risk increasing)',
    'LOG_ODDS': 'Log odds ratio (0 --> no effect; above 0 --> A1 is risk increasing)',
    'INFO': 'INFO score (imputation quality; higher --> better imputation)',
    'FRQ': 'Allele frequency',
    'SIGNED_SUMSTAT': 'Directional summary statistic as specified by --signed-sumstats.',
    'NSTUDY': 'Number of studies in which the SNP was genotyped.'
}

numeric_cols = ['P', 'N', 'N_CAS', 'N_CON', 'Z', 'OR', 'BETA', 'LOG_ODDS', 'INFO', 'FRQ', 'SIGNED_SUMSTAT', 'NSTUDY']

__version__ = '2.0.1'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* LD Score Regression (LDSC)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane\n"
MASTHEAD += "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f


class Logger(object):
    '''
    Lightweight logging.
    TODO: replace with logging module

    '''
    def __init__(self, fh):
        self.log_fh = open(fh, 'wb')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        TODO: print to file
        '''
        ##does change below still print to file?
        #print(msg, file=self.log_fh)
        print(msg)

def get_compression(fh):
    '''
    Read filename suffixes and figure out whether it is gzipped,bzip2'ed or not compressed
    '''
    if fh.endswith('gz'):
        compression = 'gzip'
        openfunc = gzip.open
    elif fh.endswith('bz2'):
        compression = 'bz2'
        openfunc = bz2.BZ2File
    else:
        openfunc = open
        compression = None

    return openfunc, compression

def read_header(fh):
    '''Read the first line of a file and returns a list with the column names.'''
    (openfunc, compression) = get_compression(fh)
    return [x.rstrip('\n') for x in openfunc(fh).readline().split()]

def parse_flag_cnames(log, args):
    '''
    Parse flags that specify how to interpret nonstandard column names.

    flag_cnames is a dict that maps (cleaned) arguments to internal column names
    '''
    cname_options = [
        [args.nstudy, 'NSTUDY', '--nstudy'],
        [args.snp, 'SNP', '--snp'],
        [args.N_col, 'N', '--N'],
        [args.N_cas_col, 'N_CAS', '--N-cas-col'],
        [args.N_con_col, 'N_CON', '--N-con-col'],
        [args.a1, 'A1', '--a1'],
        [args.a2, 'A2', '--a2'],
        [args.p, 'P', '--P'],
        [args.frq, 'FRQ', '--nstudy'],
        [args.info, 'INFO', '--info']
    ]
    flag_cnames = {clean_header(x[0]): x[1]
                   for x in cname_options if x[0] is not None}
    if args.info_list:
        try:
            flag_cnames.update(
                {clean_header(x): 'INFO' for x in args.info_list.split(',')})
        except ValueError:
            log.log(
                'The argument to --info-list should be a comma-separated list of column names.')
            raise

    null_value = None
    if args.signed_sumstats:
        try:
            cname, null_value = args.signed_sumstats.split(',')
            null_value = float(null_value)
            flag_cnames[clean_header(cname)] = 'SIGNED_SUMSTAT'
        except ValueError:
            log.log(
                'The argument to --signed-sumstats should be column header comma number.')
            raise

    return [flag_cnames, null_value]


def clean_header(header):
    '''
    For cleaning file headers.
    - convert to uppercase
    - replace dashes '-' with underscores '_'
    - replace dots '.' (as in R) with underscores '_'
    - remove newlines ('\n')
    '''
    return header.upper().replace('-', '_').replace('.', '_').replace('\n', '')

def get_cname_map(flag, default, ignore):
    '''
    Figure out which column names to use.

    Priority is
    (1) ignore everything in ignore
    (2) use everything in flags that is not in ignore
    (3) use everything in default that is not in ignore or in flags

    The keys of flag are cleaned. The entries of ignore are not cleaned. The keys of defualt
    are cleaned. But all equality is modulo clean_header().

    '''
    clean_ignore = [clean_header(x) for x in ignore]
    cname_map = {x: flag[x] for x in flag if x not in clean_ignore}
    cname_map.update(
        {x: default[x] for x in default if x not in clean_ignore + list(flag.keys())})
    return cname_map


def filter_pvals(P, log, args):
    '''Remove out-of-bounds P-values'''
    ii = (P > 0) & (P <= 1)
    bad_p = (~ii).sum()
    if bad_p > 0:
        msg = 'WARNING: {N} SNPs had P outside of (0,1]. The P column may be mislabeled.'
        log.log(msg.format(N=bad_p))

    return ii


def filter_info(info, log, args):
    '''Remove INFO < args.info_min (default 0.9) and complain about out-of-bounds INFO.'''
    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= args.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = (((info > 2.0) & info.notnull()).any(axis=1) | (
            (info < 0) & info.notnull()).any(axis=1))
        ii = (info.sum(axis=1) >= args.info_min * (len(info.columns)))
    else:
        raise ValueError('Expected pd.DataFrame or pd.Series.')

    bad_info = jj.sum()
    if bad_info > 0:
        msg = 'WARNING: {N} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled.'
        log.log(msg.format(N=bad_info))

    return ii


def filter_frq(frq, log, args):
    '''
    Filter on MAF. Remove MAF < args.maf_min and out-of-bounds MAF.
    '''
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    if bad_frq > 0:
        msg = 'WARNING: {N} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled.'
        log.log(msg.format(N=bad_frq))

    frq = np.minimum(frq, 1 - frq)
    ii = frq > args.maf_min
    return ii & ~jj

_N_CHR = 22
# complementary bases
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# bases
BASES = COMPLEMENT.keys()
# true iff strand ambiguous
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                    for x in it.product(BASES, BASES)
                    if x[0] != x[1]}
# SNPS we want to keep (pairs of alleles)
VALID_SNPS = {x for x in list(map(lambda y: ''.join(y), it.product(BASES, BASES)))
              if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {x for x in list(map(lambda y: ''.join(y), it.product(VALID_SNPS, VALID_SNPS)))
                 # strand and ref match
                 if ((x[0] == x[2]) and (x[1] == x[3])) or
                 # ref match, strand flip
                 ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                 # ref flip, strand match
                 ((x[0] == x[3]) and (x[1] == x[2])) or
                 ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip
# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {''.join(x):
                ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                # strand flip
                ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                for x in MATCH_ALLELES}


def filter_alleles(a):
    '''Remove alleles that do not describe strand-unambiguous SNPs'''
    return a.isin(VALID_SNPS)

def parse_dat(dat_gen, convert_colname, merge_alleles, log, args):
    '''Parse and filter a sumstats file chunk-wise'''
    tot_snps = 0
    dat_list = []
    msg = 'Reading sumstats from {F} into memory {N} SNPs at a time.'
    log.log(msg.format(F=args.sumstats, N=int(args.chunksize)))
    drops = {'NA': 0, 'P': 0, 'INFO': 0,
             'FRQ': 0, 'A': 0, 'SNP': 0, 'MERGE': 0}
    for block_num, dat in enumerate(dat_gen):
        sys.stdout.write('.')
        tot_snps += len(dat)
        old = len(dat)
        dat = dat.dropna(axis=0, how="any", subset=list(filter(lambda x: x != 'INFO', dat.columns))).reset_index(drop=True)
        drops['NA'] += old - len(dat)
        dat.columns = map(lambda x: convert_colname[x], dat.columns)

        wrong_types = [c for c in dat.columns if c in numeric_cols and not np.issubdtype(dat[c].dtype, np.number)]
        if len(wrong_types) > 0:
            raise ValueError('Columns {} are expected to be numeric'.format(wrong_types))

        ii = np.array([True for i in range(len(dat))])
        if args.merge_alleles:
            old = ii.sum()
            ii = dat.SNP.isin(merge_alleles.SNP)
            drops['MERGE'] += old - ii.sum()
            if ii.sum() == 0:
                continue

            dat = dat[ii].reset_index(drop=True)
            ii = np.array([True for i in range(len(dat))])

        if 'INFO' in dat.columns:
            old = ii.sum()
            ii &= filter_info(dat['INFO'], log, args)
            new = ii.sum()
            drops['INFO'] += old - new
            old = new

        if 'FRQ' in dat.columns:
            old = ii.sum()
            ii &= filter_frq(dat['FRQ'], log, args)
            new = ii.sum()
            drops['FRQ'] += old - new
            old = new

        old = ii.sum()
        if args.keep_maf:
            dat.drop(
                [x for x in ['INFO'] if x in dat.columns], inplace=True, axis=1)
        else:
            dat.drop(
                [x for x in ['INFO', 'FRQ'] if x in dat.columns], inplace=True, axis=1)
        ii &= filter_pvals(dat.P, log, args)
        new = ii.sum()
        drops['P'] += old - new
        old = new
        if not args.no_alleles:
            dat.A1 = dat.A1.str.upper()
            dat.A2 = dat.A2.str.upper()
            ii &= filter_alleles(dat.A1 + dat.A2)
            new = ii.sum()
            drops['A'] += old - new
            old = new

        if ii.sum() == 0:
            continue

        dat_list.append(dat[ii].reset_index(drop=True))

    sys.stdout.write(' done\n')
    dat = pd.concat(dat_list, axis=0).reset_index(drop=True)
    msg = 'Read {N} SNPs from --sumstats file.\n'.format(N=tot_snps)
    if args.merge_alleles:
        msg += 'Removed {N} SNPs not in --merge-alleles.\n'.format(
            N=drops['MERGE'])

    msg += 'Removed {N} SNPs with missing values.\n'.format(N=drops['NA'])
    msg += 'Removed {N} SNPs with INFO <= {I}.\n'.format(
        N=drops['INFO'], I=args.info_min)
    msg += 'Removed {N} SNPs with MAF <= {M}.\n'.format(
        N=drops['FRQ'], M=args.maf_min)
    msg += 'Removed {N} SNPs with out-of-bounds p-values.\n'.format(
        N=drops['P'])
    msg += 'Removed {N} variants that were not SNPs or were strand-ambiguous.\n'.format(
        N=drops['A'])
    msg += '{N} SNPs remain.'.format(N=len(dat))
    log.log(msg)
    return dat

def process_n(dat, args, log):
    '''Determine sample size from --N* flags or N* columns. Filter out low N SNPs.s'''
    if all(i in dat.columns for i in ['N_CAS', 'N_CON']):
        N = dat.N_CAS + dat.N_CON
        P = dat.N_CAS / N
        dat['N'] = N * P / P[N == N.max()].mean()
        dat.drop(['N_CAS', 'N_CON'], inplace=True, axis=1)
        # NB no filtering on N done here -- that is done in the next code block

    if 'N' in dat.columns:
        n_min = args.n_min if args.n_min else dat.N.quantile(0.9) / 1.5
        old = len(dat)
        dat = dat[dat.N >= n_min].reset_index(drop=True)
        new = len(dat)
        log.log('Removed {M} SNPs with N < {MIN} ({N} SNPs remain).'.format(
            M=old - new, N=new, MIN=n_min))

    elif 'NSTUDY' in dat.columns and 'N' not in dat.columns:
        nstudy_min = args.nstudy_min if args.nstudy_min else dat.NSTUDY.max()
        old = len(dat)
        dat = dat[dat.NSTUDY >= nstudy_min].drop(
            ['NSTUDY'], axis=1).reset_index(drop=True)
        new = len(dat)
        log.log('Removed {M} SNPs with NSTUDY < {MIN} ({N} SNPs remain).'.format(
            M=old - new, N=new, MIN=nstudy_min))

    if 'N' not in dat.columns:
        if args.N:
            dat['N'] = args.N
            log.log('Using N = {N}'.format(N=args.N))
        elif args.N_cas and args.N_con:
            dat['N'] = args.N_cas + args.N_con
            if args.daner is None:
                msg = 'Using N_cas = {N1}; N_con = {N2}'
                log.log(msg.format(N1=args.N_cas, N2=args.N_con))
        else:
            raise ValueError('Cannot determine N. This message indicates a bug.\n'
                             'N should have been checked earlier in the program.')

    return dat


def p_to_z(P, N):
    '''Convert P-value and N to standardized beta.'''
    return np.sqrt(chi2.isf(P, 1))

def check_median(x, expected_median, tolerance, name):
    '''Check that median(x) is within tolerance of expected_median.'''
    m = np.median(x)
    if np.abs(m - expected_median) > tolerance:
        msg = 'WARNING: median value of {F} is {V} (should be close to {M}). This column may be mislabeled.'
        raise ValueError(msg.format(F=name, M=expected_median, V=round(m, 2)))
    else:
        msg = 'Median value of {F} was {C}, which seems sensible.'.format(
            C=m, F=name)

    return msg

def allele_merge(dat, alleles, log):
    '''
    WARNING: dat now contains a bunch of NA's~
    Note: dat now has the same SNPs in the same order as --merge alleles.
    '''
    dat = pd.merge(
        alleles, dat, how='left', on='SNP', sort=False).reset_index(drop=True)
    ii = dat.A1.notnull()
    a1234 = dat.A1[ii] + dat.A2[ii] + dat.MA[ii]
    match = a1234.apply(lambda y: y in MATCH_ALLELES)
    jj = pd.Series(np.zeros(len(dat), dtype=bool))
    jj[ii] = match
    old = ii.sum()
    n_mismatch = (~match).sum()
    if n_mismatch < old:
        log.log('Removed {M} SNPs whose alleles did not match --merge-alleles ({N} SNPs remain).'.format(M=n_mismatch,
                                                                                                         N=old - n_mismatch))
    else:
        raise ValueError(
            'All SNPs have alleles that do not match --merge-alleles.')

    dat.loc[~jj.astype('bool'), [i for i in dat.columns if i != 'SNP']] = float('nan')
    dat.drop(['MA'], axis=1, inplace=True)
    return dat




######################


args = parser.parse_args('')
dir(args)
print(args)
args.merge_alleles = 'C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\w_hm3.snplist'
args.sumstats = 'C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\iq.sumstats\\iq_df_22.txt'
#args.sumstats = iq_df_22
args.out = 'iq_22'
p=True
##cant have both B and Z in file
args.ignore='B'
args.a1 = 'EffectAllele'
args.a2 = 'OtherAllele'

#run object setting lined at beginning of munge_sumstats.

# set p = False for testing in order to prevent printing
START_TIME = time.time()
log = Logger(args.out + '.log')
try:
    if p:
        defaults = vars(parser.parse_args(''))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += "Call: \n"
        header += './munge_sumstats.py \\\n'
        options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)

    file_cnames = read_header(args.sumstats)  # note keys not cleaned
    print(file_cnames)
    flag_cnames, signed_sumstat_null = parse_flag_cnames(log, args)
    print(flag_cnames)

    ### saved as list doesn't work
    ignore_cnames = [clean_header(x) for x in args.ignore.split(',')]    # remove LOG_ODDS, BETA, Z, OR from the default list
    mod_default_cnames = default_cnames
    print(mod_default_cnames)    
    #2 dictionaries and a list
    cname_map = get_cname_map(flag_cnames, mod_default_cnames, ignore_cnames)

    cname_translation = {x: cname_map[clean_header(x)] for x in file_cnames if
                         clean_header(x) in cname_map}  # note keys not cleaned
    print(cname_translation)
    cname_description = {
        x: describe_cname[cname_translation[x]] for x in cname_translation}
    print(cname_translation)    
    #print(args.signed_sumstats)
    #print(args.a1_inc)
    #print(args)
    
    sign_cnames = [
        x for x in cname_translation if cname_translation[x] in null_values]
    print(sign_cnames)
    if len(sign_cnames) > 1:
        raise ValueError(
            'Too many signed sumstat columns. Specify which to ignore with the --ignore flag.')
    if len(sign_cnames) == 0:
        raise ValueError(
            'Could not find a signed summary statistic column.')
    sign_cname = sign_cnames[0]
    signed_sumstat_null = null_values[cname_translation[sign_cname]]
    cname_translation[sign_cname] = 'SIGNED_SUMSTAT'    
    

   # print(len(sign_cnames))

    #print(sign_cnames)
    #print(sign_cname)
    #print(signed_sumstat_null)
    #print(cname_translation)

    # check that we have all the columns we need
    if not args.a1_inc:
        req_cols = ['SNP', 'P', 'SIGNED_SUMSTAT']
    else:
        req_cols = ['SNP', 'P']

    #print(req_cols)
    #print(cname_translation.values())
    #print(cname_translation)
    for c in req_cols:
        if c not in cname_translation.values():
            raise ValueError('Could not find {C} column.'.format(C=c))

    # check aren't any duplicated column names in mapping
    for field in cname_translation:
        numk = file_cnames.count(field)
        if numk > 1:
            raise ValueError('Found {num} columns named {C}'.format(C=field,num=str(numk)))

    # check multiple different column names don't map to same data field
    
    for head in cname_translation.values():
        numc = list(cname_translation.values()).count(head)
        if numc > 1:
            raise ValueError('Found {num} different {C} columns'.format(C=head,num=str(numc)))

    if (not args.N) and (not (args.N_cas and args.N_con)) and ('N' not in cname_translation.values()) and\
            (any(x not in cname_translation.values() for x in ['N_CAS', 'N_CON'])):
        raise ValueError('Could not determine N.')
    if ('N' in cname_translation.values() or all(x in cname_translation.values() for x in ['N_CAS', 'N_CON']))\
            and 'NSTUDY' in cname_translation.values():
        nstudy = [
            x for x in cname_translation if cname_translation[x] == 'NSTUDY']
        for x in nstudy:
            del cname_translation[x]
    if not args.no_alleles and not all(x in cname_translation.values() for x in ['A1', 'A2']):
        raise ValueError('Could not find A1/A2 columns.')

    log.log('Interpreting column names as follows:')
    log.log('\n'.join([x + ':\t' + cname_description[x]
                       for x in cname_description]) + '\n')

    if args.merge_alleles:
        log.log(
            'Reading list of SNPs for allele merge from {F}'.format(F=args.merge_alleles))
        (openfunc, compression) = get_compression(args.merge_alleles)
        merge_alleles = pd.read_csv(args.merge_alleles, compression=compression, header=0,
                                    sep='\s+', na_values='.')
        if any(x not in merge_alleles.columns for x in ["SNP", "A1", "A2"]):
            raise ValueError(
                '--merge-alleles must have columns SNP, A1, A2.')

        log.log(
            'Read {N} SNPs for allele merge.'.format(N=len(merge_alleles)))
        merge_alleles['MA'] = (
            merge_alleles.A1 + merge_alleles.A2).apply(lambda y: y.upper())
        merge_alleles.drop(
            [x for x in merge_alleles.columns if x not in ['SNP', 'MA']], axis=1, inplace=True)
    else:
        merge_alleles = None
    
    #print(merge_alleles.head())
    (openfunc, compression) = get_compression(args.sumstats)
    #print(openfunc)
    #print(compression)

    # figure out which columns are going to involve sign information, so we can ensure
    # they're read as floats
    signed_sumstat_cols = [k for k,v in cname_translation.items() if v=='SIGNED_SUMSTAT']
    #print(signed_sumstat_cols)
    dat_gen = pd.read_csv(args.sumstats, sep='\s+', header=0,
            compression=compression, usecols=cname_translation.keys(),
            na_values=['.', 'NA'], iterator=True, chunksize=args.chunksize,
            dtype={c:np.float64 for c in signed_sumstat_cols})
    #print(dat_gen)
    
    dat = parse_dat(dat_gen, cname_translation, merge_alleles, log, args)
    if len(dat) == 0:
        raise ValueError('After applying filters, no SNPs remain.')
    
    #print(dat.tail(500))    
    #print(dat.Z.min())
    #print(dat.Z.min())

    old = len(dat)
    #print(old)
    dat = dat.drop_duplicates(subset='SNP').reset_index(drop=True)
    new = len(dat)
    log.log('Removed {M} SNPs with duplicated rs numbers ({N} SNPs remain).'.format(
        M=old - new, N=new))
    # filtering on N cannot be done chunkwise
    dat = process_n(dat, args, log)
    dat.P = p_to_z(dat.P, dat.N)
    dat.rename(columns={'P': 'Z'}, inplace=True)
    if not args.a1_inc:
        log.log(
            check_median(dat.SIGNED_SUMSTAT, signed_sumstat_null, 0.1, sign_cname))
        dat.Z *= (-1) ** (dat.SIGNED_SUMSTAT < signed_sumstat_null)
        dat.drop('SIGNED_SUMSTAT', inplace=True, axis=1)
    # do this last so we don't have to worry about NA values in the rest of
    # the program
    '''
    FutureWarning: Setting an item of incompatible dtype is deprecated and will raise an error in a future version of pandas. 
    Value '[True True True ... True True True]' has dtype incompatible with bool, please explicitly cast to a compatible dtype first.
    '''
    if args.merge_alleles:
        dat = allele_merge(dat, merge_alleles, log)

    out_fname = args.out + '.sumstats'
    print_colnames = [
        c for c in dat.columns if c in ['SNP', 'N', 'Z', 'A1', 'A2']]
    if args.keep_maf and 'FRQ' in dat.columns:
        print_colnames.append('FRQ')
    msg = 'Writing summary statistics for {M} SNPs ({N} with nonmissing beta) to {F}.'
    log.log(
        msg.format(M=len(dat), F=out_fname + '.gz', N=dat.N.notnull().sum()))
    if p:
        dat.to_csv(out_fname + '.gz', sep="\t", index=False,
                   columns=print_colnames, float_format='%.3f', compression = 'gzip')

    log.log('\nMetadata:')
    CHISQ = (dat.Z ** 2)
    #print(CHISQ)
    mean_chisq = CHISQ.mean()
    log.log('Mean chi^2 = ' + str(round(mean_chisq, 3)))
    #print(mean_chisq)
    if mean_chisq < 1.02:
        log.log("WARNING: mean chi^2 may be too small.")

    log.log('Lambda GC = ' + str(round(CHISQ.median() / 0.4549, 3)))
    log.log('Max chi^2 = ' + str(round(CHISQ.max(), 3)))
    log.log('{N} Genome-wide significant SNPs (some may have been removed by filtering).'.format(N=(CHISQ
                                                                                                    > 29).sum()))
    #return dat

except Exception:
    log.log('\nERROR converting summary statistics:\n')
    ex_type, ex, tb = sys.exc_info()
    log.log(traceback.format_exc(ex))
    raise
finally:
    log.log('\nConversion finished at {T}'.format(T=time.ctime()))
    log.log('Total time elapsed: {T}'.format(
        T=sec_to_str(round(time.time() - START_TIME, 2))))
    
    
#dat.to_csv('test.tsv', sep='\t')    

###update allele_merge

###