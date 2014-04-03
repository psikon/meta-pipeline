#!/usr/bin/env python

#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# IMPORTS
import sys, os
from argparse import ArgumentParser
import shlex
import subprocess

# GLOBAL VARIABLES
trimmomatic = 'java -Xmx8G -jar /home/psehnert/software/Trimmomatic/trimmomatic-0.32.jar'
    
def main(argv = None):

    # Setup argument parser
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

    if __name__ == '__main__':

        try:
            os.makedirs(args.output)
        except OSError:
            # if dir exists and is dir go ahead
            if not os.path.isdir(args.output):
                raise


        readname1 = args.input[0].split('/')[-1].split('.')[0]
        readname2 = args.input[1].split('/')[-1].split('.')[0]
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
        trimming_summary = [int(s) for s in trim.stderr.read().split() if s.isdigit()]
        sys.stdout.write('Input Reads: %d          \n\
                          Both Surviving: %d (%d\%) \n\
                          Forward only: %d (%d\%)   \n\
                          Reverse only: %d (%d\%)   \n\
                          Filtered out: %d (%d\%)   \n' % 
                        (trimming_summary[-5],
                        trimming_summary[-4], round(trimming_summary[-4]*100/trimming_summary[-5],2),
                        trimming_summary[-3], round(trimming_summary[-3]*100/trimming_summary[-5],2),
                        trimming_summary[-2], round(trimming_summary[-2]*100/trimming_summary[-5],2),
                        trimming_summary[-1], 1))
        # Input Read Pairs: 1400909 Both Surviving: 1072089 (76,53%) Forward Only Surviving: 178601 (12,75%) Reverse Only Surviving: 38126 (2,72%) Dropped: 112093 (8,00%)
        # TrimmomaticPE: Completed successfully

        input = [args.output + os.sep + readname1 + '.trimmed.fastq', args.output + os.sep + readname2 + '.trimmed.fastq']
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
        filter_summary = [int(s) for s in filter.stderr.read().split() if s.isdigit()]
        sys.stdout.write('Input Reads: %d          \n\
                          Both Surviving: %d (%d%) \n\
                          Forward only: %d (%d%)   \n\
                          Reverse only: %d (%d%)   \n\
                          Filtered out: %d (%d%)   \n' % 
                        (filter_summary[-5], round(filter_summary[-4]*100/filter_summary[-5],2),
                        filter_summary[-4], round(filter_summary[-3]*100/filter_summary[-5],2),
                        filter_summary[-3], round(filter_summary[-2]*100/filter_summary[-5],2),
                        filter_summary[-2], round(filter_summary[-1]*100/filter_summary[-5],2),
                        filter_summary[-1], 1))
sys.exit(main())