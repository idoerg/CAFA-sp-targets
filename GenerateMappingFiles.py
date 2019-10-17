#!/usr/bin/env python
import glob
import re
import os
import sys

"""
Maps UniprotKB/Swissprot IDs from the cafa target files to the CAFA IDs
Files to be mapped should have the format <somestring>.taxid.(fasta|tfa)
"""

def make_mapping_files(targetdir="TargetFiles",mapdir="MappingFiles"):
    filepat = re.compile(r'%s\/[^\.]+\.([0-9]+)\.(tfa|fasta)' % targetdir)
    if not os.path.exists(mapdir):
        os.makedirs(mapdir) 
    for infilename in glob.glob(targetdir+"/*"):
        mid = re.match(filepat, infilename)
        if not mid:
            continue
        outfilename = "%s/mapping.%s.map" % (mapdir,str(mid.group(1)))
        fh1 = open(infilename, "r")  
        fh2 = open(outfilename, "w")
        for line in fh1:
            if line[0] == ">":
                fh2.write(str("\t".join(line[1:].split()))+"\n")
        fh2.close()
        fh1.close()
if __name__ == "__main__":
    make_mapping_files(sys.argv[1], sys.argv[2])
