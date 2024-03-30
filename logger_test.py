# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:09:29 2024

@author: antsh
"""
import time, argparse #sys, traceback,
from functools import reduce

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
        self.log_fh = open(fh, 'a',encoding='utf-8')

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.

        '''
        ##does change below still print to file?
        print(msg, file=self.log_fh)
        print(msg)

class Logger2(object):
    def __init__(self, fh):
        self.log_fh = open(fh, 'a', encoding='utf-8')  # Open in text mode and specify encoding

    def log(self, msg):
        print(msg, file=self.log_fh)
        print(msg)

    def close(self):
        self.log_fh.close()  # Close the file when done
        
class Logger3(object):
    def __init__(self, fh):
        self.file_handle = fh

    def __enter__(self):
        self.log_fh = open(self.file_handle, 'a', encoding='utf-8')  # Open the file
        return self

    def log(self, msg):
        print(msg, file=self.log_fh)
        print(msg)

    def __exit__(self, exc_type, exc_value, traceback):
        self.log_fh.close()  # Close the file when exiting the 'with' block
        

 # Any additional logging code here

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

parser = argparse.ArgumentParser()
parser.add_argument('--out', default='ldsc', type=str,
    help='Output filename prefix. If --out is not set, LDSC will use ldsc as the '
    'defualt output filename prefix.')
args = parser.parse_args('')
args.out="C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\03output\\log_test_3"


log = Logger(args.out+'.log')
log2 = Logger2(args.out+'.log2')
log3 = Logger3(args.out+'.log3')

header = MASTHEAD
log.log(header)
log2.log(header)
log3.log(header)
with Logger3(args.out+'.log3') as logger:
 logger.log(header)

 
log.log_fh
#header += "Call: \n"
#header += './ldsc.py \\\n'
#options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
#header += '\n'.join(options).replace('True','').replace('False','')
header = header[0:-1]+'\n'
log.log(header)
log.log('Beginning analysis at {T}'.format(T=time.ctime()))
start_time = time.time()

log.log('Analysis finished at {T}'.format(T=time.ctime()) )
time_elapsed = round(time.time()-start_time,2)
log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))



class Logger(object):
    def __init__(self, fh):
        self.log_fh = open(fh, 'w')

    def log(self, msg):
        print(msg, file=self.log_fh)
        print(msg)

out="C:\\Users\\antsh\\Documents\\sullivan\\ldsc_exome\\03output\\log_test_2"
log = Logger(out+'.log')

header = "some text"
log.log(header)