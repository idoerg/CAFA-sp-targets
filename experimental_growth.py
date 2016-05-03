#!/usr/bin/env python
#
# Growth of experimentally annotated proteins in SwissProt
# in different species / ontologies
#
import sys
from Bio import SwissProt as sp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

from filter_sp_species import go_ec_filter

exp_ec = ['EXP','IDA','IPI','IMP','IGI','IEP']

def n_go_ec(sp_handle, taxon_id, allowed_ec=exp_ec):
    # For each GO namespace, number of proteins annotated with allowed
    # terms. Uses experimental evidence code, by default.
    # 
    ncount = {}
    ncount['BPO'] = 0
    ncount['MFO'] = 0
    ncount['CCO'] = 0
    for inrec in sp.parse(sp_handle):
        if taxon_id in inrec.taxonomy_id:
            in_allowed = go_ec_filter(inrec,allowed_ec=exp_ec)
            for onto in in_allowed:
                if in_allowed[onto]:
                    ncount[onto] += 1
    return ncount

def count_exp_prots(sp_file_list, taxon_id):
    ncount_sp = {}
    for sp_file in sp_file_list:
        sp_handle = open(sp_file)
        ncount = n_go_ec(sp_handle, taxon_id)
        ncount_sp[sp_file] = ncount
    return ncount_sp

def print_count_exp_prots(sp_file_list, taxon_id):
    ncount_sp = count_exp_prots(sp_file_list,taxon_id)
    print("SP_version\tMFO\tBPO\tCCO")
    for sp_file in sp_file_list:
        mfo_count = ncount_sp[sp_file]['MFO']
        bpo_count = ncount_sp[sp_file]['BPO']
        cco_count = ncount_sp[sp_file]['CCO']
        print("%s\t%d\t%d\t%d" % (sp_file,
                    mfo_count, bpo_count, cco_count))
def read_sp_file_list(infile):
    sp_file_list = []
    with open(infile) as f:
        for inline in f:
            sp_file_list.append(inline.strip())
    return sp_file_list

if __name__ == '__main__':
    sp_file_list = read_sp_file_list(sys.argv[1])
    print_count_exp_prots(sp_file_list, sys.argv[2])
