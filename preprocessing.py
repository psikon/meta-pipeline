#!/usr/bin/env python

'''
litte wrapper script for automated quality improvement of Illumina fastq files in 3 steps:
- quality based trimming of reads
- length filtering
- dereplication/removing of duplicates
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
    parser.add_argument('input', nargs = '+', action = 'store', 
                        help = 'single or paired input files in <fastq> format')

    
    args = parser.parse_args()

    # import path of executable scripts from settings.conf
    conf = ConfigParser.SafeConfigParser()
    conf.read('settings.conf')
    trimmomatic = 'java -Xmx12G -jar ' + conf.get('Executable','trimmomatic')

    if __name__ == '__main__':
        # create output dir
        try:
            os.makedirs(args.output)
        except OSError:
            # if dir exists and is dir go ahead
            if not os.path.isdir(args.output):
                raise
        # check path of executable
        try:
          os.path.exists(conf.get('Executable','trimmomatic'))
        except:
          sys.stderr.write("Executable not found! Please check settings.conf")
          sys.exit(1)

        # create names for output files
        readname1 = args.input[0].split('/')[-1].split('.')[0]
        readname2 = args.input[1].split('/')[-1].split('.')[0]

#################################################################################################
        
        sys.stdout.write('Starting quality based trimming with args:\n\
                         LEADING: %d\n\
                         TRAILING: %d\n\
                         SLIDING_WINDOW: %s\n' % (args.leading, args.trailing, args.sliding_window))
        # quality based trimming of 3' and 5' ends with sliding window algoriythm
        trim = subprocess.Popen(shlex.split('%s PE -threads %d -phred33 -trimlog %s %s %s %s %s %s LEADING:%d TRAILING:%d SLIDINGWINDOW:%s' % 
                                            (trimmomatic,
                                            args.threads,
                                            args.output + os.sep + readname1 + '.trim.log',
                                            ' '.join(str(i)for i in args.input),
                                            args.output + os.sep + readname1 + '.trimmed.fastq',
                                            args.output + os.sep + readname1 + '.unpaired.fastq',
                                            args.output + os.sep + readname2 + '.trimmed.fastq',
                                            args.output + os.sep + readname2 + '.unpaired.fastq',
                                            args.leading, 
                                            args.trailing, 
                                            args.sliding_window)),
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
        # create new input files for next step
        input = [args.output + os.sep + readname1 + '.trimmed.fastq', args.output + os.sep + readname2 + '.trimmed.fastq']
        
###################################################################################################################################
        
        sys.stdout.write('Starting quality based trimming with args:\nMINLEN: %d\n' % (args.minlength))
        # length filtering of reads
        filter = subprocess.Popen(shlex.split('%s PE -threads %d -phred33 -trimlog %s %s %s %s %s %s MINLEN:%d' %
                                             (trimmomatic,
                                             args.threads,
                                             args.output + os.sep + readname1 + '.filtered.log',
                                             ' '.join(str(i)for i in input),
                                             args.output + os.sep + readname1 + '.filtered.fastq',
                                             args.output + os.sep + readname1 + '.short.fastq',
                                             args.output + os.sep + readname2 + '.filtered.fastq',
                                             args.output + os.sep + readname2 + '.short.fastq',
                                             args.minlength)),
                                            stderr = subprocess.PIPE)
        filter.wait()
        # parse cmd output
        filter_summary = [int(s) for s in filter.stderr.read().split() if s.isdigit()]
        # new cmd output
        sys.stdout.write('Input Reads: %d          \n\
                          Both Surviving: %d - %5.2f%%  \n\
                          Forward only: %d - %5.2f%%    \n\
                          Reverse only: %d  - %5.2f%%   \n\
                          Filtered out: %d - %5.2f%%    \n' % (filter_summary[-5], 
                                                    filter_summary[-4], 
                                                    0.0 if filter_summary[-4] == 0 else round(filter_summary[-4]*100/filter_summary[-5],2),
                                                    filter_summary[-3], 
                                                    0.0 if filter_summary[-3] == 0 else round(filter_summary[-3]*100/filter_summary[-5],2),
                                                    filter_summary[-2], 
                                                    0.0 if filter_summary[-2] == 0 else round(filter_summary[-2]*100/filter_summary[-5],2),
                                                    filter_summary[-1], 
                                                    0.0 if filter_summary[-1] == 0 else round(filter_summary[-1]*100/filter_summary[-5],2))
                          )
sys.exit(main())