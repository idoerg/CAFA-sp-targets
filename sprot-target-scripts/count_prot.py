#!/usr/bin/env python
import sys
import argparse
from collections import defaultdict
from Bio.UniProtGOA import UniProtGOA as upg
if __name__ == '__main__':
    protdict = defaultdict(int)
    for infile in sys.argv[1:]:
        s = 0
        for i in upg.gafbyproteiniterator(open(infile)):
            s += 1    
        print infile, s
