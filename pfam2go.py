#!/usr/bin/env python
'''
script to create out of a pfam annotation file a new file with go annotation
'''
#@author: Philipp Sehnert
#@contact: philipp.sehnert[a]gmail.com

# IMPORTS
import sys, os
from argparse import ArgumentParser
from itertools import islice

# GLOBAL VARIABLES

def create_go_table(mapping_file):
    # init the go dictionary
    go_table = {}
    sys.stdout.write('Import Pfam2GO Index from %s\n' % (mapping_file))
    # open go annotation file
    with open(mapping_file,'r') as f:
        # remove first line with header
        next(f)
        next(f)
        # iterate over all lines and create a key:value pair for every line in the file
        for line in f:
            # remove \n
            line = line.strip()
            # split line into fields
            fields = line.split('\t')
            # create key: value pair: key = PfamID:value = rest of line 
            go_table[fields[0]] = fields[1:]
        f.close()
    sys.stdout.write('Loaded %d GO annotations.\n' % (len(go_table)))
    return(go_table)

def read_pfam_file(pfam_file):
    # init the pfam dictionary
    pfam_table = {}
    sys.stdout.write('Load pfam file: %s\n' % (pfam_file))
    # open the pfam annotation file
    with open(pfam_file,'r') as f:
        # iterate over the file line by line and create key:value pairs
        for line in f:
            # ignore commented lines at start
            if not line.startswith('#'):
                # remove \n
                line = line.strip()
                # split the line into fields
                fields = line.split()
                # key = PfamID
                key = fields[3].split('.',1)[0]
                # create key:value pair with value = rest of line
                pfam_table[key] = fields[0:]
        f.close()
    sys.stdout.write('Successfully imported %d pfam annotations\n' % (len(pfam_table)))
    return(pfam_table)

def compare_keys(pfam_table, go_table):
    sys.stdout.write('Annotate Pfam with GOs.\n')
    # init count and array for matching key combinations
    count = 0
    key_list = []
    # iterate over pfam annotation keys
    for key in pfam_table.keys():
        # compare actual key with key from go table
        if key in go_table.keys():
            # raise count if combination was found and add key to array
            count +=1
            key_list.append(key)

    sys.stdout.write('Successfully annotated %d/%d pfams\n' % (count, len(pfam_table)))
    return(key_list)

def create_output_table(key_list, pfam_table, go_table, outputfile):
    # open an output file and create lines with following fields
    # pfam target_seq GO GO name GO origin e-value score 
    with open(outputfile,'w') as f:
        # write header
        f.write('#pfam \t target_seq \t GO-ID \t GO-Desc \t GO-Tree \t e-value, \t score\n')
        # for every key in the array create a new line and write it to file
        for key in key_list:
            pfam = pfam_table.get(key)
            go = go_table.get(key)
            line = [key,pfam[0],go[0], go[1], go[2], pfam[4], pfam[5]]
            # create tab seperated line in file
            f.write('\t'.join(map(str,line)))
            f.write('\n')

def main(argv = None):

    # Setup cmd interface
    parser = ArgumentParser(description = '%s -- map GO terms to pfam annotation' % 
                            (os.path.basename(sys.argv[0])),
                            epilog = 'created by Philipp Sehnert')
    parser.add_argument('-i', dest = 'pfam', required = True,
                        help = 'location of pfam input file')
    parser.add_argument('-o', dest = 'output', required = True, 
                        help = 'destination for output')
    parser.add_argument('-m', dest='mapping', required = True,
                        help = ' loaction of pfam2go mapping file')
    # parse arguments from cmd interface
    args = parser.parse_args()

    if __name__ == '__main__':
        go = create_go_table(args.mapping)
        pfam = read_pfam_file(args.pfam)
        key_list = compare_keys(pfam,go)
        create_output_table(key_list,pfam,go,args.output)

sys.exit(main())