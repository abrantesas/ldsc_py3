# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 07:55:01 2024
Partitioned Heritability
@author: antsh
"""

import ldscore.ldscore as ld
import ldscore.parse as ps
#import ldscore.sumstats as sumstats
import ldscore.regressions as reg
import numpy as np
import pandas as pd
from subprocess import call
from itertools import product
import time, sys, traceback, argparse
from functools import reduce
import gzip
import copy
import itertools as it
import os
from scipy import stats

try:
    x = pd.DataFrame({'A': [1, 2, 3]})
    x.sort_values(by='A')
except AttributeError:
    raise ImportError('LDSC requires pandas version >= 0.17.0')

__version__ = '2.0.1'
MASTHEAD = "*********************************************************************\n"
MASTHEAD += "* LD Score Regression (LDSC)\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2014-2019 Brendan Bulik-Sullivan and Hilary Finucane\n"
MASTHEAD += "* Broad Institute of MIT and Harvard / MIT Department of Mathematics\n"
MASTHEAD += "* Performance enhancement and python3 adaption by Anthony Abrantes 2024\n"
MASTHEAD += "* University of North Carolina-Chapel Hill\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "*********************************************************************\n"
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.precision', 4)
pd.set_option('max_colwidth',1000)
np.set_printoptions(linewidth=1000)
np.set_printoptions(precision=4)


parser = argparse.ArgumentParser()
parser.add_argument('--out', default='ldsc', type=str,
    help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
    'defualt output filename prefix.')
# Basic LD Score Estimation Flags'
parser.add_argument('--bfile', default=None, type=str,
    help='Prefix for Plink .bed/.bim/.fam file')
parser.add_argument('--l2', default=False, action='store_true',
    help='Estimate l2. Compatible with both jackknife and non-jackknife.')
# Filtering / Data Management for LD Score
parser.add_argument('--extract', default=None, type=str,
    help='File with SNPs to include in LD Score estimation. '
    'The file should contain one SNP ID per row.')
parser.add_argument('--keep', default=None, type=str,
    help='File with individuals to include in LD Score estimation. '
    'The file should contain one individual ID per row.')
parser.add_argument('--ld-wind-snps', default=None, type=int,
    help='Specify the window size to be used for estimating LD Scores in units of '
    '# of SNPs. You can only specify one --ld-wind-* option.')
parser.add_argument('--ld-wind-kb', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'kilobase-pairs (kb). You can only specify one --ld-wind-* option.')
parser.add_argument('--ld-wind-cm', default=None, type=float,
    help='Specify the window size to be used for estimating LD Scores in units of '
    'centiMorgans (cM). You can only specify one --ld-wind-* option.')
parser.add_argument('--print-snps', default=None, type=str,
    help='This flag tells LDSC to only print LD Scores for the SNPs listed '
    '(one ID per row) in PRINT_SNPS. The sum r^2 will still include SNPs not in '
    'PRINT_SNPs. This is useful for reducing the number of LD Scores that have to be '
    'read into memory when estimating h2 or rg.' )
# Fancy LD Score Estimation Flags
parser.add_argument('--annot', default=None, type=str,
    help='Filename prefix for annotation file for partitioned LD Score estimation. '
    'LDSC will automatically append .annot or .annot.gz to the filename prefix. '
    'See docs/file_formats_ld for a definition of the .annot format.')
parser.add_argument('--thin-annot', action='store_true', default=False,
    help='This flag says your annot files have only annotations, with no SNP, CM, CHR, BP columns.')
parser.add_argument('--cts-bin', default=None, type=str,
    help='This flag tells LDSC to compute partitioned LD Scores, where the partition '
    'is defined by cutting one or several continuous variable[s] into bins. '
    'The argument to this flag should be the name of a single file or a comma-separated '
    'list of files. The file format is two columns, with SNP IDs in the first column '
    'and the continuous variable in the second column. ')
parser.add_argument('--cts-breaks', default=None, type=str,
    help='Use this flag to specify names for the continuous variables cut into bins '
    'with --cts-bin. For each continuous variable, specify breaks as a comma-separated '
    'list of breakpoints, and separate the breakpoints for each variable with an x. '
    'For example, if binning on MAF and distance to gene (in kb), '
    'you might set --cts-breaks 0.1,0.25,0.4x10,100,1000 ')
parser.add_argument('--cts-names', default=None, type=str,
    help='Use this flag to specify names for the continuous variables cut into bins '
    'with --cts-bin. The argument to this flag should be a comma-separated list of '
    'names. For example, if binning on DAF and distance to gene, you might set '
    '--cts-bin DAF,DIST_TO_GENE ')
parser.add_argument('--per-allele', default=False, action='store_true',
    help='Setting this flag causes LDSC to compute per-allele LD Scores, '
    'i.e., \ell_j := \sum_k p_k(1-p_k)r^2_{jk}, where p_k denotes the MAF '
    'of SNP j. ')
parser.add_argument('--pq-exp', default=None, type=float,
    help='Setting this flag causes LDSC to compute LD Scores with the given scale factor, '
    'i.e., \ell_j := \sum_k (p_k(1-p_k))^a r^2_{jk}, where p_k denotes the MAF '
    'of SNP j and a is the argument to --pq-exp. ')
parser.add_argument('--no-print-annot', default=False, action='store_true',
    help='By defualt, seting --cts-bin or --cts-bin-add causes LDSC to print '
    'the resulting annot matrix. Setting --no-print-annot tells LDSC not '
    'to print the annot matrix. ')
parser.add_argument('--maf', default=None, type=float,
    help='Minor allele frequency lower bound. Default is MAF > 0.')
# Basic Flags for Working with Variance Components
parser.add_argument('--h2', default=None, type=str,
    help='Filename for a .sumstats[.gz] file for one-phenotype LD Score regression. '
    '--h2 requires at minimum also setting the --ref-ld and --w-ld flags.')
parser.add_argument('--h2-cts', default=None, type=str,
    help='Filename for a .sumstats[.gz] file for cell-type-specific analysis. '
    '--h2-cts requires the --ref-ld-chr, --w-ld, and --ref-ld-chr-cts flags.')
parser.add_argument('--rg', default=None, type=str,
    help='Comma-separated list of prefixes of .chisq filed for genetic correlation estimation.')
parser.add_argument('--ref-ld', default=None, type=str,
    help='Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD '
    'Score regression. '
    'LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix.')
parser.add_argument('--ref-ld-chr', default=None, type=str,
    help='Same as --ref-ld, but will automatically concatenate .l2.ldscore files split '
    'across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz '
    'to the filename prefix. If the filename prefix contains the symbol @, LDSC will '
    'replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome '
    'numbers to the end of the filename prefix.'
    'Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz'
    'Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz')
parser.add_argument('--w-ld', default=None, type=str,
    help='Filename prefix for file with LD Scores with sum r^2 taken over SNPs included '
    'in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz.')
parser.add_argument('--w-ld-chr', default=None, type=str,
    help='Same as --w-ld, but will read files split into 22 chromosomes in the same '
    'manner as --ref-ld-chr.')
parser.add_argument('--overlap-annot', default=False, action='store_true',
    help='This flag informs LDSC that the partitioned LD Scores were generates using an '
    'annot matrix with overlapping categories (i.e., not all row sums equal 1), '
    'and prevents LDSC from displaying output that is meaningless with overlapping categories.')
parser.add_argument('--print-coefficients',default=False,action='store_true',
    help='when categories are overlapping, print coefficients as well as heritabilities.')
parser.add_argument('--frqfile', type=str,
    help='For use with --overlap-annot. Provides allele frequencies to prune to common '
    'snps if --not-M-5-50 is not set.')
parser.add_argument('--frqfile-chr', type=str,
    help='Prefix for --frqfile files split over chromosome.')
parser.add_argument('--no-intercept', action='store_true',
    help = 'If used with --h2, this constrains the LD Score regression intercept to equal '
    '1. If used with --rg, this constrains the LD Score regression intercepts for the h2 '
    'estimates to be one and the intercept for the genetic covariance estimate to be zero.')
parser.add_argument('--intercept-h2', action='store', default=None,
    help = 'Intercepts for constrained-intercept single-trait LD Score regression.')
parser.add_argument('--intercept-gencov', action='store', default=None,
    help = 'Intercepts for constrained-intercept cross-trait LD Score regression.'
    ' Must have same length as --rg. The first entry is ignored.')
parser.add_argument('--M', default=None, type=str,
    help='# of SNPs (if you don\'t want to use the .l2.M files that came with your .l2.ldscore.gz files)')
parser.add_argument('--two-step', default=None, type=float,
    help='Test statistic bound for use with the two-step estimator. Not compatible with --no-intercept and --constrain-intercept.')
parser.add_argument('--chisq-max', default=None, type=float,
    help='Max chi^2.')
parser.add_argument('--ref-ld-chr-cts', default=None, type=str,
    help='Name of a file that has a list of file name prefixes for cell-type-specific analysis.')
parser.add_argument('--print-all-cts', action='store_true', default=False)

# Flags for both LD Score estimation and h2/gencor estimation
parser.add_argument('--print-cov', default=False, action='store_true',
    help='For use with --h2/--rg. This flag tells LDSC to print the '
    'covaraince matrix of the estimates.')
parser.add_argument('--print-delete-vals', default=False, action='store_true',
    help='If this flag is set, LDSC will print the block jackknife delete-values ('
    'i.e., the regression coefficeints estimated from the data with a block removed). '
    'The delete-values are formatted as a matrix with (# of jackknife blocks) rows and '
    '(# of LD Scores) columns.')
# Flags you should almost never use
parser.add_argument('--chunk-size', default=50, type=int,
    help='Chunk size for LD Score calculation. Use the default.')
parser.add_argument('--pickle', default=False, action='store_true',
    help='Store .l2.ldscore files as pickles instead of gzipped tab-delimited text.')
parser.add_argument('--yes-really', default=False, action='store_true',
    help='Yes, I really want to compute whole-chromosome LD Score.')
parser.add_argument('--invert-anyway', default=False, action='store_true',
    help="Force LDSC to attempt to invert ill-conditioned matrices.")
parser.add_argument('--n-blocks', default=200, type=int,
    help='Number of block jackknife blocks.')
parser.add_argument('--not-M-5-50', default=False, action='store_true',
    help='This flag tells LDSC to use the .l2.M file instead of the .l2.M_5_50 file.')
parser.add_argument('--return-silly-things', default=False, action='store_true',
    help='Force ldsc to return silly genetic correlation estimates.')
parser.add_argument('--no-check-alleles', default=False, action='store_true',
    help='For rg estimation, skip checking whether the alleles match. This check is '
    'redundant for pairs of chisq files generated using munge_sumstats.py and the '
    'same argument to the --merge-alleles flag.')
# transform to liability scale
parser.add_argument('--samp-prev',default=None,
    help='Sample prevalence of binary phenotype (for conversion to liability scale).')
parser.add_argument('--pop-prev',default=None,
    help='Population prevalence of binary phenotype (for conversion to liability scale).')
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


def _remove_dtype(x):
    '''Removes dtype: float64 and dtype: int64 from pandas printouts'''
    x = str(x)
    x = x.replace('\ndtype: int64', '')
    x = x.replace('\ndtype: float64', '')
    return x


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

        '''
        ##does change below still print to file?
        #print(msg, file=self.log_fh)
        print(msg)


def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = ps.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        print(f(c, len(x.IDList)))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            print(f(c, len_merged_list))
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list

def annot_sort_key(s):
    '''For use with --cts-bin. Fixes weird pandas crosstab column order.'''
    if type(s) == tuple:
        s = [x.split('_')[0] for x in s]
        ## added list to py3 changed functionality of map
        s = list(map(lambda x: float(x) if x != 'min' else -float('inf'), s))
    else:  # type(s) = str:
        s = s.split('_')[0]
        if s == 'min':
            s = float('-inf')
        else:
            s = float(s)

    return s

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


def _splitp(fstr):
    flist = fstr.split(',')
    flist = [os.path.expanduser(os.path.expandvars(x)) for x in flist]
    return flist


def _select_and_log(x, ii, log, msg):
    '''Fiter down to rows that are True in ii. Log # of SNPs removed.'''
    new_len = ii.sum()
    if new_len == 0:
        raise ValueError(msg.format(N=0))
    else:
        x = x[ii]
        log.log(msg.format(N=new_len))
    return x


def smart_merge(x, y):
    '''Check if SNP columns are equal. If so, save time by using concat instead of merge.'''
    if len(x) == len(y) and (x.index == y.index).all() and (x.SNP == y.SNP).all():
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True).drop('SNP', 1)
        out = pd.concat([x, y], axis=1)
    else:
        out = pd.merge(x, y, how='inner', on='SNP')
    return out


def _read_ref_ld(args, log):
    '''Read reference LD Scores.'''
    ref_ld = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                   'reference panel LD Score', ps.ldscore_fromlist)
    log.log(
        'Read reference panel LD Scores for {N} SNPs.'.format(N=len(ref_ld)))
    return ref_ld


def _read_annot(args, log):
    '''Read annot matrix.'''
    try:
        if args.ref_ld is not None:
            overlap_matrix, M_tot = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                                          'annot matrix', ps.annot, frqfile=args.frqfile)
        elif args.ref_ld_chr is not None:
            overlap_matrix, M_tot = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                                      'annot matrix', ps.annot, frqfile=args.frqfile_chr)
    except Exception:
        log.log('Error parsing .annot file.')
        raise

    return overlap_matrix, M_tot


def _read_M(args, log, n_annot):
    '''Read M (--M, --M-file, etc).'''
    if args.M:
        try:
            M_annot = [float(x) for x in _splitp(args.M)]
        except ValueError as e:
            raise ValueError('Could not cast --M to float: ' + str(e.args))
    else:
        if args.ref_ld:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld_chr), _N_CHR, common=(not args.not_M_5_50))

    try:
        M_annot = np.array(M_annot).reshape((1, n_annot))
    except ValueError as e:
        raise ValueError(
            '# terms in --M must match # of LD Scores in --ref-ld.\n' + str(e.args))

    return M_annot


def _read_w_ld(args, log):
    '''Read regression SNP LD.'''
    if (args.w_ld and ',' in args.w_ld) or (args.w_ld_chr and ',' in args.w_ld_chr):
        raise ValueError(
            '--w-ld must point to a single fileset (no commas allowed).')
    w_ld = _read_chr_split_files(args.w_ld_chr, args.w_ld, log,
                                 'regression weight LD Score', ps.ldscore_fromlist)
    if len(w_ld.columns) != 2:
        raise ValueError('--w-ld may only have one LD Score column.')
    w_ld.columns = ['SNP', 'LD_weights']  # prevent colname conflicts w/ ref ld
    log.log(
        'Read regression weight LD Scores for {N} SNPs.'.format(N=len(w_ld)))
    return w_ld


def _read_chr_split_files(chr_arg, not_chr_arg, log, noun, parsefunc, **kwargs):
    '''Read files split across 22 chromosomes (annot, ref_ld, w_ld).'''
    try:
        if not_chr_arg:
            log.log('Reading {N} from {F} ... ({p})'.format(N=noun, F=not_chr_arg, p=parsefunc.__name__))
            out = parsefunc(_splitp(not_chr_arg), **kwargs)
        elif chr_arg:
            f = ps.sub_chr(chr_arg, '[1-22]')
            log.log('Reading {N} from {F} ... ({p})'.format(N=noun, F=f, p=parsefunc.__name__))
            out = parsefunc(_splitp(chr_arg), _N_CHR, **kwargs)
    except ValueError as e:
        log.log('Error parsing {N}.'.format(N=noun))
        raise e

    return out


def _read_sumstats(args, log, fh, alleles=False, dropna=False):
    '''Parse summary statistics.'''
    log.log('Reading summary statistics from {S} ...'.format(S=fh))
    sumstats = ps.sumstats(fh, alleles=alleles, dropna=dropna)
    log_msg = 'Read summary statistics for {N} SNPs.'
    log.log(log_msg.format(N=len(sumstats)))
    m = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset='SNP')
    if m > len(sumstats):
        log.log(
            'Dropped {M} SNPs with duplicated rs numbers.'.format(M=m - len(sumstats)))

    return sumstats


def _check_ld_condnum(args, log, ref_ld):
    '''Check condition number of LD Score matrix.'''
    if len(ref_ld.shape) >= 2:
        cond_num = int(np.linalg.cond(ref_ld))
        if cond_num > 100000:
            if args.invert_anyway:
                warn = "WARNING: LD Score matrix condition number is {C}. "
                warn += "Inverting anyway because the --invert-anyway flag is set."
                log.log(warn.format(C=cond_num))
            else:
                warn = "WARNING: LD Score matrix condition number is {C}. "
                warn += "Remove collinear LD Scores. "
                raise ValueError(warn.format(C=cond_num))


def _check_variance(log, M_annot, ref_ld):
    '''Remove zero-variance LD Scores.'''
    ii = ref_ld.iloc[:, 1:].var() == 0  # NB there is a SNP column here
    if ii.all():
        raise ValueError('All LD Scores have zero variance.')
    else:
        log.log('Removing partitioned LD Scores with zero variance.')
        ii_snp = np.array([True] + list(~ii))
        ii_m = np.array(~ii)
        ref_ld = ref_ld.iloc[:, ii_snp]
        M_annot = M_annot[:, ii_m]

    return M_annot, ref_ld, ii


def _warn_length(log, sumstats):
    if len(sumstats) < 200000:
        log.log(
            'WARNING: number of SNPs less than 200k; this is almost always bad.')


def _print_cov(ldscore_reg, ofh, log):
    '''Prints covariance matrix of slopes.'''
    log.log(
        'Printing covariance matrix of the estimates to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.coef_cov)


def _print_delete_values(ldscore_reg, ofh, log):
    '''Prints block jackknife delete-k values'''
    log.log('Printing block jackknife delete values to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.tot_delete_values)

def _print_part_delete_values(ldscore_reg, ofh, log):
    '''Prints partitioned block jackknife delete-k values'''
    log.log('Printing partitioned block jackknife delete values to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.part_delete_values)


def _merge_and_log(ld, sumstats, noun, log):
    '''Wrap smart merge with log messages about # of SNPs.'''
    sumstats = smart_merge(ld, sumstats)
    msg = 'After merging with {F}, {N} SNPs remain.'
    if len(sumstats) == 0:
        raise ValueError(msg.format(N=len(sumstats), F=noun))
    else:
        log.log(msg.format(N=len(sumstats), F=noun))

    return sumstats


def _read_ld_sumstats(args, log, fh, alleles=False, dropna=True):
    sumstats = _read_sumstats(args, log, fh, alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args, log)
    n_annot = len(ref_ld.columns) - 1
    M_annot = _read_M(args, log, n_annot)
    M_annot, ref_ld, novar_cols = _check_variance(log, M_annot, ref_ld)
    w_ld = _read_w_ld(args, log)
    sumstats = _merge_and_log(ref_ld, sumstats, 'reference panel LD', log)
    sumstats = _merge_and_log(sumstats, w_ld, 'regression SNP LD', log)
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.columns[1:len(ref_ld.columns)]
    return M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols

def cell_type_specific(args, log):
    '''Cell type specific analysis'''
    args = copy.deepcopy(args)
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1

    M_annot_all_regr, w_ld_cname, ref_ld_cnames_all_regr, sumstats, novar_cols = \
            _read_ld_sumstats(args, log, args.h2_cts)
    #M_tot is assigned but not used anywhere?
    #M_tot = np.sum(M_annot_all_regr)
    _check_ld_condnum(args, log, ref_ld_cnames_all_regr)
    _warn_length(log, sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    if args.chisq_max is None:
        chisq_max = max(0.001*sumstats.N.max(), 80)
    else:
        chisq_max = args.chisq_max

    ii = np.ravel(sumstats.Z**2 < chisq_max)
    sumstats = sumstats.iloc[ii, :]
    log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
            C=chisq_max, N=np.sum(ii), M=n_snp-np.sum(ii)))
    n_snp = np.sum(ii)  # lambdas are late-binding, so this works
    ref_ld_all_regr = np.array(sumstats[ref_ld_cnames_all_regr]).reshape((len(sumstats),-1))
    chisq = np.array(sumstats.Z**2)
    keep_snps = sumstats[['SNP']]

    s = lambda x: np.array(x).reshape((n_snp, 1))
    results_columns = ['Name', 'Coefficient', 'Coefficient_std_error', 'Coefficient_P_value']
    results_data = []
    for (name, ct_ld_chr) in [x.split() for x in open(args.ref_ld_chr_cts).readlines()]:
        ref_ld_cts_allsnps = _read_chr_split_files(ct_ld_chr, None, log,
                                   'cts reference panel LD Score', ps.ldscore_fromlist)
        log.log('Performing regression.')
        ref_ld_cts = np.array(pd.merge(keep_snps, ref_ld_cts_allsnps, on='SNP', how='left').iloc[:,1:])
        if np.any(np.isnan(ref_ld_cts)):
            raise ValueError ('Missing some LD scores from cts files. Are you sure all SNPs in ref-ld-chr are also in ref-ld-chr-cts')

        ref_ld = np.hstack([ref_ld_cts, ref_ld_all_regr])
        M_cts = ps.M_fromlist(
                _splitp(ct_ld_chr), _N_CHR, common=(not args.not_M_5_50))
        M_annot = np.hstack([M_cts, M_annot_all_regr])
        #what is this?
        hsqhat = reg.Hsq(s(chisq), ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N),
                     M_annot, n_blocks=n_blocks, intercept=args.intercept_h2,
                     twostep=None, old_weights=True)
        coef, coef_se = hsqhat.coef[0], hsqhat.coef_se[0]
        results_data.append((name, coef, coef_se, stats.norm.sf(coef/coef_se)))
        if args.print_all_cts:
            for i in range(1, len(ct_ld_chr.split(','))):
                coef, coef_se = hsqhat.coef[i], hsqhat.coef_se[i]
                results_data.append((name+'_'+str(i), coef, coef_se, stats.norm.sf(coef/coef_se)))


    df_results = pd.DataFrame(data = results_data, columns = results_columns)
    df_results.sort_values(by = 'Coefficient_P_value', inplace=True)
    df_results.to_csv(args.out+'.cell_type_results.txt', sep='\t', index=False)
    log.log('Results printed to '+args.out+'.cell_type_results.txt')


def estimate_rg(args, log):
    '''Estimate rg between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)
    rg_paths, rg_files = _parse_rg(args.rg)
    n_pheno = len(rg_paths)
    f = lambda x: _split_or_none(x, n_pheno)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev = list(map(f,(args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev)))
    list(map(lambda x: _check_arg_len(x, n_pheno), ((args.intercept_h2, '--intercept-h2'),
                                               (args.intercept_gencov, '--intercept-gencov'),
                                               (args.samp_prev, '--samp-prev'),
                                               (args.pop_prev, '--pop-prev'))))
    if args.no_intercept:
        args.intercept_h2 = [1 for _ in range(n_pheno)]
        args.intercept_gencov = [0 for _ in range(n_pheno)]
    p1 = rg_paths[0]
    out_prefix = args.out + rg_files[0]
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, p1,
                                                                        alleles=True, dropna=True)
    RG = []
    n_annot = M_annot.shape[1]
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    if args.two_step is not None:
        log.log('Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))

    for i, p2 in enumerate(rg_paths[1:n_pheno]):
        log.log(
            'Computing rg for phenotype {I}/{N}'.format(I=i + 2, N=len(rg_paths)))
        try:
            loop = _read_other_sumstats(args, log, p2, sumstats, ref_ld_cnames)
            rghat = _rg(loop, args, log, M_annot, ref_ld_cnames, w_ld_cname, i)
            RG.append(rghat)
            _print_gencor(args, log, rghat, ref_ld_cnames, i, rg_paths, i == 0)
            out_prefix_loop = out_prefix + '_' + rg_files[i + 1]
            if args.print_cov:
                _print_rg_cov(rghat, out_prefix_loop, log)
            if args.print_delete_vals:
                _print_rg_delete_values(rghat, out_prefix_loop, log)

        except Exception:  # keep going if phenotype 50/100 causes an error
            msg = 'ERROR computing rg for phenotype {I}/{N}, from file {F}.'
            log.log(msg.format(I=i + 2, N=len(rg_paths), F=rg_paths[i + 1]))
            ex_type, ex, tb = sys.exc_info()
            log.log(traceback.format_exc(ex) + '\n')
            if len(RG) <= i:  # if exception raised before appending to RG
                RG.append(None)

    log.log('\nSummary of Genetic Correlation Results\n' +
            _get_rg_table(rg_paths, RG, args))
    return RG


def _read_other_sumstats(args, log, p2, sumstats, ref_ld_cnames):
    loop = _read_sumstats(args, log, p2, alleles=True, dropna=False)
    loop = _merge_sumstats_sumstats(args, sumstats, loop, log)
    loop = loop.dropna(how='any')
    alleles = loop.A1 + loop.A2 + loop.A1x + loop.A2x
    if not args.no_check_alleles:
        loop = _select_and_log(loop, _filter_alleles(alleles), log,
                               '{N} SNPs with valid alleles.')
        loop['Z2'] = _align_alleles(loop.Z2, alleles)

    loop = loop.drop(['A1', 'A1x', 'A2', 'A2x'], axis=1)
    _check_ld_condnum(args, log, loop[ref_ld_cnames])
    _warn_length(log, loop)
    return loop


def _get_rg_table(rg_paths, RG, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, 'NA')
    x = pd.DataFrame()
    x['p1'] = [rg_paths[0] for i in range(1, len(rg_paths))]
    x['p2'] = rg_paths[1:len(rg_paths)]
    x['rg'] = list(map(t('rg_ratio'), RG))
    x['se'] = list(map(t('rg_se'), RG))
    x['z'] = list(map(t('z'), RG))
    x['p'] = list(map(t('p'), RG))
    #typo below? i or it?
    if args.samp_prev is not None and \
            args.pop_prev is not None and \
            all((i is not None for i in args.samp_prev)) and \
            all((it is not None for it in args.pop_prev)):

        c = list(map(lambda x, y: reg.h2_obs_to_liab(1, x, y), args.samp_prev[1:], args.pop_prev[1:]))
        x['h2_liab'] = list(map(lambda x, y: x * y, c, list(map(t('tot')), list(map(t('hsq2')), RG))))
        x['h2_liab_se'] = list(map(lambda x, y: x * y, c, list(map(t('tot_se')), list(map(t('hsq2')), RG))))
    else:
        x['h2_obs'] = list(map(t('tot'), list(map(t('hsq2'), RG))))
        x['h2_obs_se'] = list(map(t('tot_se'), list(map(t('hsq2'), RG))))

    x['h2_int'] = list(map(t('intercept'), list(map(t('hsq2'), RG))))
    x['h2_int_se'] = list(map(t('intercept_se'), list(map(t('hsq2'), RG))))
    x['gcov_int'] = list(map(t('intercept'), list(map(t('gencov'), RG))))
    x['gcov_int_se'] = list(map(t('intercept_se'), list(map(t('gencov'), RG))))
    return x.to_string(header=True, index=False) + '\n'


def _print_gencor(args, log, rghat, ref_ld_cnames, i, rg_paths, print_hsq1):
    l = lambda x: x + ''.join(['-' for i in range(len(x.replace('\n', '')))])
    P = [args.samp_prev[0], args.samp_prev[i + 1]]
    K = [args.pop_prev[0], args.pop_prev[i + 1]]
    if args.samp_prev is None and args.pop_prev is None:
        args.samp_prev = [None, None]
        args.pop_prev = [None, None]
    if print_hsq1:
        log.log(l('\nHeritability of phenotype 1\n'))
        log.log(rghat.hsq1.summary(ref_ld_cnames, P=P[0], K=K[0]))

    log.log(
        l('\nHeritability of phenotype {I}/{N}\n'.format(I=i + 2, N=len(rg_paths))))
    log.log(rghat.hsq2.summary(ref_ld_cnames, P=P[1], K=K[1]))
    log.log(l('\nGenetic Covariance\n'))
    log.log(rghat.gencov.summary(ref_ld_cnames, P=P, K=K))
    log.log(l('\nGenetic Correlation\n'))
    log.log(rghat.summary() + '\n')


def _merge_sumstats_sumstats(args, sumstats1, sumstats2, log):
    '''Merge two sets of summary statistics.'''
    sumstats1.rename(columns={'N': 'N1', 'Z': 'Z1'}, inplace=True)
    sumstats2.rename(
        columns={'A1': 'A1x', 'A2': 'A2x', 'N': 'N2', 'Z': 'Z2'}, inplace=True)
    x = _merge_and_log(sumstats1, sumstats2, 'summary statistics', log)
    return x


def _filter_alleles(alleles):
    '''Remove bad variants (mismatched alleles, non-SNPs, strand ambiguous).'''
    ii = alleles.apply(lambda y: y in MATCH_ALLELES)
    #ii = alleles.apply(lambda y: y in MATCH_ALLELES)
    return ii


def _align_alleles(z, alleles):
    '''Align Z1 and Z2 to same choice of ref allele (allowing for strand flip).'''
    try:
        z *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
        #z *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
    except KeyError as e:
        msg = 'Incompatible alleles in .sumstats files: %s. ' % e.args
        msg += 'Did you forget to use --merge-alleles with munge_sumstats.py?'
        raise KeyError(msg)
    return z


def _rg(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname, i):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    s = lambda x: np.array(x).reshape((n_snp, 1))
    if args.chisq_max is not None:
        ii = sumstats.Z1**2*sumstats.Z2**2 < args.chisq_max**2
        n_snp = np.sum(ii)  # lambdas are late binding, so this works
        sumstats = sumstats[ii]
    n_blocks = min(args.n_blocks, n_snp)
    #removed deprecated as_matrix
    ref_ld = sumstats[ref_ld_cnames].values
    intercepts = [args.intercept_h2[0], args.intercept_h2[
        i + 1], args.intercept_gencov[i + 1]]
    rghat = reg.RG(s(sumstats.Z1), s(sumstats.Z2),
                   ref_ld, s(sumstats[w_ld_cname]), s(
                       sumstats.N1), s(sumstats.N2), M_annot,
                   intercept_hsq1=intercepts[0], intercept_hsq2=intercepts[1],
                   intercept_gencov=intercepts[2], n_blocks=n_blocks, twostep=args.two_step)

    return rghat


def _parse_rg(rg):
    '''Parse args.rg.'''
    rg_paths = _splitp(rg)
    rg_files = [x.split('/')[-1] for x in rg_paths]
    if len(rg_paths) < 2:
        raise ValueError(
            'Must specify at least two phenotypes for rg estimation.')

    return rg_paths, rg_files


def _print_rg_delete_values(rg, fh, log):
    '''Print block jackknife delete values.'''
    _print_delete_values(rg.hsq1, fh + '.hsq1.delete', log)
    _print_delete_values(rg.hsq2, fh + '.hsq2.delete', log)
    _print_delete_values(rg.gencov, fh + '.gencov.delete', log)


def _print_rg_cov(rghat, fh, log):
    '''Print covariance matrix of estimates.'''
    _print_cov(rghat.hsq1, fh + '.hsq1.cov', log)
    _print_cov(rghat.hsq2, fh + '.hsq2.cov', log)
    _print_cov(rghat.gencov, fh + '.gencov.cov', log)


def _split_or_none(x, n):
    if x is not None:
        y = list(map(float, x.replace('N', '-').split(',')))
    else:
        y = [None for _ in range(n)]
    return y


def _check_arg_len(x, n):
    x, m = x
    if len(x) != n:
        raise ValueError(
            '{M} must have the same number of arguments as --rg/--h2.'.format(M=m))





annot = pd.read_csv('C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\baseline.22.annot',sep='\t')
print(annot.head())
print(annot.columns)
keepcols=['CHR','BP','SNP','CM','base','Coding_UCSC.bed','Conserved_LindbladToh.bed','CTCF_Hoffman.bed','DGF_ENCODE.bed','DHS_Trynka.bed','Enhancer_Hoffman.bed','Intron_UCSC.bed','PromoterFlanking_Hoffman.bed']
annot_sub = annot[keepcols]
print(annot_sub.head())
annot_sub.to_csv('C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\baseline_sub.22.annot',sep='\t')

args = parser.parse_args('')
args.h2 = "C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\iq.sumstats\\iq_22"
args.ref_ld ="C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\03output\\baseline_sub.22"
args.w_ld = ####weights.\
args.overlap_annot = True
args.frqfile = ####1000G.mac5eur.\
args.out = "C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\03output\\Iq_annot_22"

'''Estimate h2 and partitioned h2.'''
args = copy.deepcopy(args)
if args.samp_prev is not None and args.pop_prev is not None:
    args.samp_prev, args.pop_prev = list(map(float, [args.samp_prev, args.pop_prev]))
if args.intercept_h2 is not None:
    args.intercept_h2 = float(args.intercept_h2)
if args.no_intercept:
    args.intercept_h2 = 1
M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
    args, log, args.h2)
ref_ld = np.array(sumstats[ref_ld_cnames])
_check_ld_condnum(args, log, ref_ld_cnames)
_warn_length(log, sumstats)
n_snp = len(sumstats)
n_blocks = min(n_snp, args.n_blocks)
n_annot = len(ref_ld_cnames)
chisq_max = args.chisq_max
old_weights = False
if n_annot == 1:
    if args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
else:
    old_weights = True
    if args.chisq_max is None:
        chisq_max = max(0.001*sumstats.N.max(), 80)

s = lambda x: np.array(x).reshape((n_snp, 1))
chisq = s(sumstats.Z**2)
if chisq_max is not None:
    ii = np.ravel(chisq < chisq_max)
    sumstats = sumstats.iloc[ii, :]
    log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
            C=chisq_max, N=np.sum(ii), M=n_snp-np.sum(ii)))
    n_snp = np.sum(ii)  # lambdas are late-binding, so this works
    ref_ld = np.array(sumstats[ref_ld_cnames])
    chisq = chisq[ii].reshape((n_snp, 1))

if args.two_step is not None:
    log.log('Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))

hsqhat = reg.Hsq(chisq, ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N),
                 M_annot, n_blocks=n_blocks, intercept=args.intercept_h2,
                 twostep=args.two_step, old_weights=old_weights)

if args.print_cov:
    _print_cov(hsqhat, args.out + '.cov', log)
if args.print_delete_vals:
    _print_delete_values(hsqhat, args.out + '.delete', log)
    _print_part_delete_values(hsqhat, args.out + '.part_delete', log)

log.log(hsqhat.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, overlap = args.overlap_annot))
if args.overlap_annot:
    overlap_matrix, M_tot = _read_annot(args, log)

    # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
    df_results = hsqhat._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
    df_results.to_csv(args.out+'.results', sep="\t", index=False)
    log.log('Results printed to '+args.out+'.results')

return hsqhat

