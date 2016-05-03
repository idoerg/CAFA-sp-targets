#!/usr/bin/env python

import sys
import glob
targets_written = {}
def cafa_sp_table(inhandle,outhandle):
    for inrec in inhandle:
        if inrec[0]=='>':
            cafa_id, sp_id = inrec[1:].split()
            if not (cafa_id in targets_written):
                outhandle.write("%s" % inrec[1:] )
                targets_written[cafa_id] = None
    return
 
def all_cafa_sp_table(midfix,inpaths):
    outhandle = open("cafa2sp-%s.txt" % midfix,"w")
    for inpath in inpaths:
        cafa_sp_table(open(inpath),outhandle)

    outhandle.close()

if __name__ == "__main__":
    midfix = sys.argv[1]
    inpaths = glob.glob(sys.argv[2])
    all_cafa_sp_table(midfix,inpaths)

