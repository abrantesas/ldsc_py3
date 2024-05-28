# ldsc_py3 LDSC (LD SCore) for python 3

## Making changes to ldsc so it will run on py3 and optimize for speed. 

### Currently repicates results for:
ldscores  \
heritability  \
genetic correlation  \
partitioned ldscores  \
partitioned heritability

### What is next:
timing and profiling \
ldscore computation spends a lot of time in numpy functions \

TODO: fix path sys.path.append(r' C:/Users/antsh/Documents/sullivan/ldsc_exome/git/ldsc_py3/ldsc_py3')

TODO: allele FRQ file for chr1 \
TODO: profile in python2, python3, and python3 with Numba \

speed up ldscore computation  \
numba or something like emeraLD? \

test and get celltype analysis to work  \
verify munge_sumstats.py works  \

### What is different:
python3 syntax  \
fixed some typos? and possible bugs  \
removed requirement of gzip installed locally  \
implemented py3 native logging \
created version that runs without bedtools or system call to gzip (i.e. can now be used on a Windows system) \
Python 3 Windows beta version available \
 
### Data

sumstats are in /nas/depts/007/sullilab/shared/gwas_sumstats/updated_sumstats/scz2022_eur/analysis/ldsc \
and also locally in ldsc_exome/01data/scz_2022_eur \
FILE: ldsc.sumstats.gz \

### Run

python ldsc.py --bfile C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\01data\\toy_ldsc\\1000G_Phase3_plinkfiles\\1000G_EUR_Phase3_plink\\1000G.EUR.QC.1 --l2 --ld-wind-cm 1 --out C:\\Users\\antsh\Documents\\sullivan\\ldsc_exome\\03output\\ldscores\\1
python -m cProfile -o /work/users/a/b/abrantes/LDSC/03analysis/test_ldsc_py2/ldscore_annot/ldscore_1.prof ldsc.py --l2 --bfile /work/users/a/b/abrantes/LDSC/01data/test_ldsc_py2/plink/1000G.EUR.QC.1 --ld-wind-cm 1 --annot /nas/depts/007/sullilab/shared/partitioned_LDSC/1000G_EUR_Phase3_baseline_annot/baseline.1.annot --out /work/users/a/b/abrantes/LDSC/03analysis/test_ldsc_py2/ldscore_annot/baseline_1 
