import glob
import re
import os

#os.mkdir("MappingFiles/")
for name in glob.glob('TargetFiles/*'):
    fh1 = open(name, "r")  
    mid = re.match(r'TargetFiles\/target\.([0-9]+)\.fasta', name)
    #print(mid.group(1))
    outfilename = "MappingFiles/mapping."+str(mid.group(1))+".map"
    fh2 = open(outfilename, "w")
    for line in fh1:
    	if line[0] == ">":
    		fh2.write(str("\t".join(line[1:].split(" "))))
fh2.close()
fh1.close()
