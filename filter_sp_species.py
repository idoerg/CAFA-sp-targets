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

def species_filter_old(sp_handle, taxon_id):
    target_id = int(taxon_id+"0000001")
    outseq_list = []
    for inrec in sp.parse(sp_handle):
        if taxon_id in inrec.taxonomy_id:
            outseq = SeqRecord(Seq(inrec.sequence),
                   id="T"+str(target_id),
                   description = "%s" %
                   (inrec.entry_name))
            outseq_list.append(outseq)
            target_id += 1
    outhandle = open("sp_species.%s.tfa" % taxon_id,"w")
    SeqIO.write(outseq_list,outhandle,"fasta")
    outhandle.close()

go_ontology_lut = {'P': 'BPO', 'C':'CCO', 'F':'MFO'}
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

allowed_ec = ['IEA','NR', 'ND', 'IC', 'NAS','TAS', 'ISS', 'ISO', 'ISA',
              'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']
def go_ec_filter(sp_rec, allowed_ec=allowed_ec,use_has_go=True):
    in_allowed_ontologies = {'MFO': True, 'BPO':True, 'CCO':True}
    if use_has_go:
        if not has_go(sp_rec):
            # print "not has go",sp_rec.entry_name
            return in_allowed_ontologies
    go_list = parse_go(sp_rec)
    for go_rec in go_list:
        if go_rec['evidence'] not in allowed_ec:
            in_allowed_ontologies[go_rec['ontology']] = False
    return in_allowed_ontologies

    
def species_filter(sp_handle, taxon_id):
    target_id = int(taxon_id+"0000001")
    outseq_list = []
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

def species_go_filter(sp_handle, taxon_id):
    target_id = int(taxon_id+"0000001")
    outhandle = {}
    outhandle['MFO'] = open("sp_species.%s.MFO.noexp.tfa" % taxon_id,"w")
    outhandle['BPO'] = open("sp_species.%s.BPO.noexp.tfa" % taxon_id,"w")
    outhandle['CCO'] = open("sp_species.%s.CCO.noexp.tfa" % taxon_id,"w")
    outhandle_all = open("sp_species.%s.all.noexp.tfa" % taxon_id,"w")
    for inrec in sp.parse(sp_handle):
        if taxon_id in inrec.taxonomy_id:
            in_allowed = go_ec_filter(inrec)
            allcount = 0
            for onto in in_allowed:
                if in_allowed[onto]:
                    allcount += 1
                    outseq = [SeqRecord(Seq(inrec.sequence),
                        id="T"+str(target_id),
                        description = "%s" %
                        (inrec.entry_name))]
                    SeqIO.write(outseq,outhandle[onto],"fasta")
            target_id += 1
            if allcount == 3:
                SeqIO.write(outseq,outhandle_all,"fasta")
    for i in outhandle:
        outhandle[i].close()
if __name__ == '__main__':
    species_filter(open(sys.argv[1]), sys.argv[2])
    species_go_filter(open(sys.argv[1]), sys.argv[2])

