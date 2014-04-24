#!/usr/bin/env python

'''
little wrapper to run flash on the paired end reads to enlength the sequences and add remaining single end
reads to the output. Also removes duplicates
ATTENTION: Remove duplicates destroys fastq and fasta header
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# global imports
import sys, os
from argparse import ArgumentParser
import subprocess
import shlex
import fileinput

# executables
FLASH = 'ext/flash'
COLLAPSER = 'ext/fastx_collapser'

def concatenation(input, outputdir, threads):
  '''wrapper for concatination of paired end reads with flash'''
  sys.stdout.write('Concatination of paired end reads ...\n')
  # Call flash on paired end reads
  concat = subprocess.Popen(shlex.split('%s -m 10 -M 200 --interleaved-output -o concat -d %s -t %d %s %s' % (FLASH,
                                                                                                              outputdir,
                                                                                                              threads,
                                                                                                              input[0],
                                                                                                              input[1])
                                      ), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  concat.wait()
  # parse all integer from output
  concat_summary = [int(s) for s in concat.stdout.read().split() if s.isdigit()]
  # create outputmsg
  msg = 'Concatiniation complete:\n\
         Input reads: %d\n\
         Concatinated reads %d \n\
         Not concationated: %d\n\
         Percentage: %d\n' % (concat_summary[-3], concat_summary[-2], concat_summary[-1], 
                              0.0 if concat_summary[-2] == 0 else round(concat_summary[-2]*100/concat_summary[-3], 4))
  # create log file and write output 
  log = outputdir + os.sep + 'flash.log'
  with open(log,'w') as log:
    log.write(msg)
  log.close()
  # print piped output on stdout
  sys.stdout.write(msg)
  # return concatinated and not concatinated files
  return [outputdir + os.sep + 'concat.extendedFrags.fastq',
          outputdir + os.sep + 'concat.notCombined.fastq']

def remove_duplicates(input, outputdir):
  '''wrapper function to call fastx_collapser on combined single reads, result will be in fasta format'''
  sys.stdout.write('Remove duplicated reads ...\n')
  # create outputs
  output = outputdir + os.sep + 'classify.nodup.fastq'
  log = outputdir + os.sep + 'no_dup.log'
  # call fastx_collapser
  duplicates = subprocess.Popen(shlex.split('%s -Q33 -v -i %s -o %s' % (COLLAPSER,
                                                                        input,
                                                                        output)),
                                stdout = subprocess.PIPE)
  duplicates.wait()
  # get piped output
  msg = duplicates.stdout.read()
  # write output to stdout ...
  sys.stdout.write(msg)
  # ... and in a logfile
  with open(log,'w') as log:
    log.write(msg)
  log.close()
  # return created file
  return output

def cat_files(input, output):
  '''combine all input file in one file'''
  sys.stdout.write('Combining all reads ...\n')
  # open output file
  with open(output, 'w') as fout:
          # append every line of every file to output file
          for line in fileinput.input(input):
            fout.write(line)
  fout.close()

  return output

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
  parser.add_argument('input', nargs = '+', action = 'store', 
                      help = 'paired end input files in <fastq> format')
  
  # get arguments of cmd
  args = parser.parse_args()
  # define first inputs
  input = args.input
  single = args.single

  if __name__ == '__main__':
    # create output dir
    try:
      os.makedirs(args.output)
    except OSError:
      # if dir exists and is dir go ahead
      if not os.path.isdir(args.output):
        raise

    try:
      # call flash
      input = concatenation(input, args.output, args.threads)
      # extend flash results with single end reads of quality control
      input.append(args.single) if args.single else None
      # combine all reads in one file
      input = cat_files(input, args.output + os.sep + 'classify.fastq')
      # remove duplicated from that file and convert to fasta
      input = remove_duplicates(input, args.output)
      sys.stdout.write('Generation of classify input complete.\nresult: %s' % (input))
    except KeyboardInterrupt:
      sys.stdout.write('\nERROR 1 : Operation cancelled by User!\n')
      sys.exit(1)
 
sys.exit(main())