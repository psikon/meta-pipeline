#!/usr/bin/env python

'''
little script to label header for qiime pick_otu.py.
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# IMPORTS
import sys, os
from argparse import ArgumentParser
import ConfigParser 
import shlex
import subprocess

# GLOBAL VARIABLES

    
def main(argv = None):

    # Setup cmd interface
    parser = ArgumentParser(description = '%s -- relabel header for pick_otu.py' % 
                            (os.path.basename(sys.argv[0])),
                            epilog = 'created by Philipp Sehnert',
    parser.add_argument('--minlength', type = int, dest = 'minlength', default = 150, required = True,
                        help = 'Drop the read if it is below a specified length')
    parser.add_argument('input', nargs = '+', action = 'store', 
                        help = 'single or paired input files in <fasta> format')

    
    args = parser.parse_args()

    if __name__ == '__main__':
      

sys.exit(main())