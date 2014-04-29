#!/usr/bin/env python

'''
litte wrapper script for automated quality improvement of Illumina fastq files in 2 steps:
- quality based trimming of reads
- length filtering
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# imports
import sys, os
from argparse import ArgumentParser
import shlex
import subprocess
import fileinput

# Executables
trimmomatic = 'java -Xmx12G -jar ext/Trimmomatic/trimmomatic-0.32.jar'

def extract_readname(item, index):
  '''extract the name of the file befor the first dot'''
  return item[index].split('/')[-1].split('.')[0]

def cat_files(input, output):
  '''combine all input file in one file'''
  # open output file
  with open(output, 'w') as fout:
          # append every line of every file to output file
          for line in fileinput.input(input):
            fout.write(line)
  fout.close()

  return output

def trimming(input, outputdir, threads, leading, trailing, sliding_window, singletons):
  '''wrapper for the trimming process with trimmomatic'''
  sys.stdout.write('Starting quality based trimming with args:\n\
                    LEADING: %d\n\
                    TRAILING: %d\n\
                    SLIDING_WINDOW: %s\n' % (leading, 
                                             trailing, 
                                             sliding_window))
  # quality based trimming of 3' and 5' ends with sliding window algorithm with trimmomatic
  trim = subprocess.Popen(shlex.split('%s PE -threads %d -phred33 -trimlog %s %s %s %s %s %s LEADING:%d TRAILING:%d SLIDINGWINDOW:%s' % 
                                      (trimmomatic,
                                       threads,
                                       outputdir + os.sep + extract_readname(input, 0) + '.trim.log',
                                       ' '.join(str(i)for i in input),
                                       outputdir + os.sep + extract_readname(input, 0) + '.trimmed.fastq',
                                       outputdir + os.sep + extract_readname(input, 0) + '.unpaired_after_trimming.fastq',
                                       outputdir + os.sep + extract_readname(input, 1) + '.trimmed.fastq',
                                       outputdir + os.sep + extract_readname(input, 1) + '.unpaired_after_trimming.fastq',
                                       leading, 
                                       trailing, 
                                       sliding_window)),
                          stderr = subprocess.PIPE)
  trim.wait()
  # parse cmd output
  trimming_summary = [int(s) for s in trim.stderr.read().split() if s.isdigit()]
  # new cmd output
  sys.stdout.write('Input Reads: %d          \n\
                    Both Surviving: %d - %5.2f%%  \n\
                    Forward only: %d - %5.2f%%    \n\
                    Reverse only: %d  - %5.2f%%   \n\
                    Filtered out: %d - %5.2f%%    \n' % (trimming_summary[-5], 
                                                         trimming_summary[-4], 
                                                         0.0 if trimming_summary[-4] == 0 else round(trimming_summary[-4]*100/trimming_summary[-5],2),
                                                         trimming_summary[-3], 
                                                         0.0 if trimming_summary[-3] == 0 else round(trimming_summary[-3]*100/trimming_summary[-5],2),
                                                         trimming_summary[-2], 
                                                         0.0 if trimming_summary[-2] == 0 else round(trimming_summary[-2]*100/trimming_summary[-5],2),
                                                         trimming_summary[-1], 
                                                         0.0 if trimming_summary[-1] == 0 else round(trimming_summary[-1]*100/trimming_summary[-5],2))
                        )
  # get successfull trimmed paired end files
  result = [outputdir + os.sep + extract_readname(input, 0) + '.trimmed.fastq', 
            outputdir + os.sep + extract_readname(input, 1) + '.trimmed.fastq']
  # get successfull trimmed but now unpaired files
  unpaired = [outputdir + os.sep + extract_readname(input, 0) + '.unpaired_after_trimming.fastq', 
              outputdir + os.sep + extract_readname(input, 1) + '.unpaired_after_trimming.fastq']
  if singletons:
    # cat forward and reverse only reads for length filtering
    single = cat_files(unpaired, 
                       outputdir + os.sep + extract_readname(input, 0) + '.single_tmp.trimmed.fastq')
  else:
    single = None
  
  # retrun paired and single results
  return [result, single]

def length_filtering_PE(input, outputdir, threads, minlength, singletons):
  '''wrapper for trimmomatic length filtering'''
  sys.stdout.write('Starting length filtering for PE with args:\nMINLEN: %d\n' % (minlength))
  # paired end length filtering of reads with trimmomatic
  filter = subprocess.Popen(shlex.split('%s PE -threads %d -phred33 -trimlog %s %s %s %s %s %s MINLEN:%d' %
                                        (trimmomatic,
                                         threads,
                                         outputdir + os.sep + extract_readname(input, 0) + '.filtered.log',
                                         ' '.join(str(i)for i in input),
                                         outputdir + os.sep + extract_readname(input, 0) + '.filtered.fastq',
                                         outputdir + os.sep + extract_readname(input, 0) + '.unpaired_after_filtering.fastq',
                                         outputdir + os.sep + extract_readname(input, 1) + '.filtered.fastq',
                                         outputdir + os.sep + extract_readname(input, 1) + '.unpaired_after_filtering.fastq',
                                         minlength)),
                            stderr = subprocess.PIPE)
  filter.wait()
  # parse cmd output
  filter_summary = [int(s) for s in filter.stderr.read().split() if s.isdigit()]
  # new cmd output
  sys.stdout.write('Input Reads: %d              \n\
                    Both Surviving: %d - %5.2f%% \n\
                    Forward only: %d - %5.2f%%   \n\
                    Reverse only: %d  - %5.2f%%  \n\
                    Filtered out: %d - %5.2f%%   \n' % (filter_summary[-5], 
                                                        filter_summary[-4], 
                                                        0.0 if filter_summary[-4] == 0 else round(filter_summary[-4]*100/filter_summary[-5],2),
                                                        filter_summary[-3], 
                                                        0.0 if filter_summary[-3] == 0 else round(filter_summary[-3]*100/filter_summary[-5],2),
                                                        filter_summary[-2], 
                                                        0.0 if filter_summary[-2] == 0 else round(filter_summary[-2]*100/filter_summary[-5],2),
                                                        filter_summary[-1], 
                                                        0.0 if filter_summary[-1] == 0 else round(filter_summary[-1]*100/filter_summary[-5],2))
                  ) 
  # get successfull filtered reads
  result = [outputdir + os.sep + extract_readname(input, 0) + '.filtered.fastq',
            outputdir + os.sep + extract_readname(input, 1) + '.filtered.fastq']
  # get unpaired reads remaining after filtering
  unpaired = [outputdir + os.sep + extract_readname(input, 0) + '.unpaired_after_filtering.fastq',
              outputdir + os.sep + extract_readname(input, 1) + '.unpaired_after_filtering.fastq']
  # if processing of singletons is switched on, then combine the unpaired reads
  if singletons:
    single = cat_files(unpaired, 
                       outputdir + os.sep + extract_readname(input, 0) + '.single_tmp.filtered.fastq')
  else:
    single = None
  # return filtered and single end reads 
  return [result, single]

def length_filtering_SE(input, outputdir, threads, minlength):
  # do length filtering for single end reads
  sys.stdout.write('Starting length filtering for SE with args:\nMINLEN: %d\n' % (minlength))
  filter = subprocess.Popen(shlex.split('%s SE -threads %d -phred33 -trimlog %s %s %s MINLEN:%d' %
                                       (trimmomatic,
                                        threads,
                                        outputdir + os.sep + extract_readname(input, 0) + '.single.log',
                                        input,
                                        outputdir + os.sep + extract_readname([input], 0) + '.single.filtered.fastq',
                                        minlength)),
                            stderr = subprocess.PIPE)
  filter.wait()
  # parse cmd output
  filter_summary = [int(s) for s in filter.stderr.read().split() if s.isdigit()]
  # new cmd output
  sys.stdout.write('Input Reads: %d\nSurviving: %d - %5.2f%%\nFiltered out: %d - %5.2f%%\n' % (filter_summary[-3], 
                                                                                               filter_summary[-2], 
                                                                                               0.0 if filter_summary[-2] == 0 else round(filter_summary[-2]*100/filter_summary[-3],2),
                                                                                               filter_summary[-1], 
                                                                                               0.0 if filter_summary[-1] == 0 else round(filter_summary[-1]*100/filter_summary[-3],2))                                                       
                  )
  # return successfull filtered single end reads
  return outputdir + os.sep + extract_readname([input], 0) + '.single.filtered.fastq'

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
  parser.add_argument('-o', dest = 'output', default = '.', required = True,
                      help = 'location for output files (default = .)')
  parser.add_argument('--leading', type = int, dest = 'leading', default = 3, required = True,
                      help = 'Cut bases off the start of a read, if below a threshold quality')
  parser.add_argument('--trailing', type = int, dest = 'trailing', default = 3, required = True,
                      help = 'Cut bases off the end of a read, if below a threshold quality')
  parser.add_argument('--sliding_window', dest = 'sliding_window', default = '4:15', required = True,
                      help = 'Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold. ')
  parser.add_argument('--minlength', type = int, dest = 'minlength', default = 150, required = True,
                      help = 'Drop the read if it is below a specified length')
  parser.add_argument('--use_no_singletons', dest = 'singletons', action = 'store_false', default = True, 
                      help = 'permit length filtering of remaining singletons reads')
  parser.add_argument('input', nargs = '+', action = 'store', 
                      help = 'single or paired input files in <fastq> format')
  # parse cmd arguments
  args = parser.parse_args()
  # define input
  input = args.input
  
  if __name__ == '__main__':
 
    # create output dir
    try:
      os.makedirs(args.output)
    except OSError:
      # if dir exists and is dir go ahead
      if not os.path.isdir(args.output):
        raise

    try:
      # start trimming process
      input = trimming(input, args.output, args.threads, 
                       args.leading, args.trailing, 
                       args.sliding_window, args.singletons)
      # seperate single end reads from trimming
      trim_single = input[1]
      # filter paired end reads for minlength
      input = length_filtering_PE(input[0], args.output, args.threads, 
                                  args.minlength, args.singletons)
      # seperate single end reads
      filtered_single = input[1]
      # combine all single end reads in one file
      all_singles_tmp = cat_files([trim_single, filtered_single], 
                               args.output + os.sep + extract_readname(input[0], 0) + '.single.fastq')
      # do a length filtereing for all remaining single end reads
      all_singles = length_filtering_SE(all_singles_tmp,
                                        args.output,
                                        args.threads,
                                        args.minlength)
       # clean up not used files
      try:
        os.remove(trim_single)
        os.remove(filtered_single)
        os.remove(all_singles_tmp)
      except:
        sys.stderr.write("Cannot cleanup completly\n")

      # give information about result files
      sys.stdout.write('Quality control complete!\nresult:\n\t%s\n\t%s\n\t%s\n' % (input[0][0],
                                                                                   input[0][1], 
                                                                                   all_singles))

    except KeyboardInterrupt:
      sys.stdout.write('\nERROR 1 : Operation cancelled by User!\n')
      sys.exit(1)

sys.exit(main())