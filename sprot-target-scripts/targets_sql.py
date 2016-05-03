
import sys
import gzip
import glob
from Bio import SwissProt as sp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def parse_go(sp_rec):
    go_list = []
    for xref in sp_rec.cross_references:
        if xref[0] == 'GO':
            go_rec = {}
            go_rec['go_id'] = xref[1]
            go_rec['ontology'] = go_ontology_lut[xref[2].split(':')[0]]
            try:
                go_rec['evidence'] = xref[3].split(':')[0]
            except IndexError:
                print sp_rec.entry_name
                print go_rec
                print xref
                raise
            go_list.append(go_rec)
    return go_list

def make_sp_table(sp_file, taxa_ids=[]):
    """
    Make a table from swissprot that includes the following fields:
    SP_ID, GO, Evidence code, Ontology
    """
    for sp_rec in sp.parse(sp_handle):
        go_list = parse_go(sp_rec)
        if not go_list:
            continue
        if taxa_ids and (not (inrec.taxonomy_id[0] in taxa_ids)):
            continue
            
        outfile.write("%s\t%s\t%s\t
