
'''
main script for meta-pipeline, that run the scripts:
	- quality_control.py
	- generate_classify.py
with given parameters.
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# global imports
import sys, os
from argparse import ArgumentParser
import subprocess
import shlex

#import generate_classify.py


def main(argv = None):

  # Setup cmd interface
  parser = ArgumentParser(description = '%s -- main script for meta-pipeline' % 
                          (os.path.basename(sys.argv[0])),
                          epilog = 'created by Philipp Sehnert',
                          add_help = True)
  parser.add_argument('--version', action = 'version', version = '%s 1.0' % 
                      (os.path.basename(sys.argv[0])))
  parser.add_argument('-t', type = int, dest = 'threads', default = 1, required = True,
                      help = 'specify the number of cpu to be used')
  parser.add_argument('-o', dest = 'output', default = '.', required = True,
  					  help = 'specify output folder')
  parser.add_argument('-s', dest = 'single',
                      help = 'include single end reads remaining after quality control')
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
    	sys.stdout.write('Running Quality Control Step\n')
    	# call quality control.py with RAW input
    	quality_control = subprocess.Popen(shlex.split('python quality_control.py -t %d -o %s --leading %d --trailing %d --sliding_window %s --minlength %s %s' % 
														(args.threads,
														 args.output + os.sep + 'quality_controled',
														 args.leading,
														 args.trailing,
														 args.sliding_window,
														 args.minlength,
														 ' '.join(input)))
    									  )
    	quality_control.wait()
      	# find quality controlled output in file structure and get usable files
      	input = []
      	# find all filtered files
      	for file in os.listdir(args.output + os.sep + 'quality_controled'):
      		if file.endswith('.filtered.fastq'):
      			input.append(file)
      	# seperate single end and paired end files
      	single = [args.output + os.sep + 'quality_controled' + os.sep + item for item in input if item.endswith('single.filtered.fastq')][0]
      	input = [args.output + os.sep + 'quality_controled' + os.sep + item for item in input if not item.endswith('single.filtered.fastq')]
      	sys.stdout.write('Running Assembly and Dereplication Step\n')
      	# call generate_classify_input.py for assembly and removing of duplicates
       	generate_classify = subprocess.Popen(shlex.split('python generate_classify_input.py -t %d -o %s -s %s %s' % 
       													(args.threads,
      		 								 			 args.output + os.sep + 'classify_input',
      		  											 single,
      		  											 ' '.join(input)))
       										)
       	generate_classify.wait()
    except KeyboardInterrupt:
      sys.stdout.write('\nERROR 1 : Operation cancelled by User!\n')
      sys.exit(1)

sys.exit(main())
