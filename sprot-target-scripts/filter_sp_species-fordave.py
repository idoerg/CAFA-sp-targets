#!/usr/bin/env python

import sys
from Bio import SwissProt as sp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
"""
    Filter SwissProt by Taxon ID

    ./filter_sp_species sp_file tax_id
"""

def has_go(sp_rec):
    retval = False
    for xref in sp_rec.cross_references:
        if xref[0] == 'GO':
            retval = True
            break
    return retval

    
def species_filter(sp_handle, taxon_id):
    target_id = int(taxon_id+"0000001")
    outhandle = open("sp_species.%s.tfa" % taxon_id,"w")
    for inrec in sp.parse(sp_handle):
        if taxon_id in inrec.taxonomy_id:
            outseq = SeqRecord(Seq(inrec.sequence),
                   id="T"+str(target_id),
                   description = "%s" %
                   (inrec.entry_name))
            outseq_list = [outseq]
            SeqIO.write(outseq_list,outhandle,"fasta")
            target_id += 1
    outhandle.close()

if __name__ == '__main__':
    species_filter(open(sys.argv[1]), sys.argv[2])

