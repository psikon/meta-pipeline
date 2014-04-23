#!/usr/bin/env python

'''
little wrapper to run flash on the paired end reads to enlength the sequences and add remaining single end
reads to the output. Also removes duplicates
ATTENTION: Remove duplicates destroys fastq and fasta header
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# IMPORTS
import sys, os
from argparse import ArgumentParser
import subprocess
import shlex

# GLOBAL VARIABLES
FLASH = 'ext/flash'
COLLAPSER = 'ext/fastx_collapser'
CONVERTER = 'ext/fastq_to_fasta'

def concatenation(input, outputdir, threads):
  concat = subprocess.Popen(shlex.split('%s -m 10 -M 200 --interleaved-output -o concat -d %s -t %d %s %s' % (FLASH,
                                                                                                        outputdir,
                                                                                                        threads,
                                                                                                        input[0],
                                                                                                        input[1])
                                      ), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  concat.wait()
  #sys.stdout.write(concat.stdout.read())
  concat_summary = [int(s) for s in concat.stdout.read().split() if s.isdigit()]
  print concat_summary


def remove_duplicates(single):
  print single

def cat_files(input,output):
  pass

def convert(input, outputdir):
  sys.stdout.write('Convert file to fasta ...')
  output = outputdir + os.sep + input.split('.')[0] + 'fasta'
  convert = subprocess.Popen(shlex('%s -i %s -o %s') % (CONVERTER,
                                                        input,
                                                        output))
  convert.wait()
  sys.stdout.write('Converted file: %s' % (output))


def main(argv = None):

  # Setup cmd interface
  parser = ArgumentParser(description = '%s -- preprocessing of paired end Illumina Reads' % 
                          (os.path.basename(sys.argv[0])),
                          epilog = 'created by Philipp Sehnert',
                          add_help = True)
  parser.add_argument('--version', action = 'version', version = '%s 1.0' % 
                      (os.path.basename(sys.argv[0])))
  parser.add_argument('-t', type = int, dest = 'threads', default = 1, required = True,
                      help = 'specify the number of cpu to be used')
  parser.add_argument('-o', dest = 'output', default = '.',
                      help = 'location for output files (default = .)')
  parser.add_argument('-s', dest = 'single',
                      help = 'include single end reads remaining after quality control')
  parser.add_argument('-f', dest = 'format', default = 'fasta', required = True,
                      choices=['fasta','fastq'], help = 'convert fastq output to fasta (for blast analysis)')
  parser.add_argument('input', nargs = '+', action = 'store', 
                      help = 'paired end input files in <fastq> format')
  
    
  args = parser.parse_args()
  input = args.input
  single = args.single

  if __name__ == '__main__':
    concatenation(input, args.output, args.threads)
    remove_duplicates(single)
    #convert()
      
sys.exit(main())