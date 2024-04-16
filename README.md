# ldsc_py3 LDSC (LD SCore) for python 3

## Making changes to ldsc so it will run on py3 and optimize for speed. 

### Currently repicates results for:
ldscores  \
heritability  \
genetic correlation  \
partitioned ldscores  \
partitioned heritability

### What is next:
timing  \
speed up ldscore computation  \
test and get celltype analysis to work  \
verify munge_sumstats.py works  \
remove reliance of bedtools and verify make_annot.py is working

### What is different:
python3 syntax  \
fixed some typos? and possible bugs  \
removed requirement of gzip installed locally  \
implemented py3 native logging