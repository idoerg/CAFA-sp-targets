#!/usr/bin/env python
import sys
import argparse
import target_prep as tp
from Bio.UniProt import GOA as upg
if __name__ == '__main__':
#    parser = argparse.ArgumentParser(description='Filter by field')
#    parser.add_argument('-o','--output')
#    parser.add_argument('-f','--field')
    outhandle = sys.stdout
    if len(sys.argv) == 5:
        outhandle = open(sys.argv[4],"w")
    outhandle.write('!gaf-version: 2.0\n')
    goodvals = {sys.argv[1]: set(sys.argv[2].split(','))}
    for inrec in upg.gafiterator(open(sys.argv[3])):
        if upg.record_has(inrec, goodvals):
            upg.writerec(inrec,outhandle)
