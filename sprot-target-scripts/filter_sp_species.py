#!/usr/bin/env python

import sys
import gzip
import glob
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

nonexp_ec = ['IEA','NR', 'ND', 'IC', 'NAS','TAS', 'ISS', 'ISO', 'ISA',
              'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']
exp_ec = ['IDA','IPI','IMP','IGI','IEP']
def go_ec_exclusive_filter(sp_rec, allowed_ec,use_has_go=True):
    """
    Filters by evidence codes. If an evidence code *is not* in the allowed
    list, then the "in_allowed_ontologies" dictionary value is set
    to False. Therefore, this is an *exclusive* filter.
    In other words, it is enough that one term in an ontology has a disallowed evidence code
    so that all the others in that ontology will be disallowed. 
    """
    in_allowed_ontologies = {'MFO': True, 'BPO':True, 'CCO':True}
    if use_has_go:
        if not has_go(sp_rec):
            # print "not has go",sp_rec.entry_name
            return in_allowed_ontologies,[]
    go_list = parse_go(sp_rec)
    ret_go_list = []
    for go_rec in go_list:
        if go_rec['evidence'] not in allowed_ec:
            in_allowed_ontologies[go_rec['ontology']] = False
        else:
            ret_go_list.append(go_rec)
    return in_allowed_ontologies, ret_go_list

def go_ec_inclusive_filter(sp_rec, allowed_ec,use_has_go=True):
    """
    Filters by evidence codes. If an evidence code is *in* the allowed
    list, then the "in_allowed_ontologies" dictionary value is set
    to True. Therefore, this is an *inclusive* filter.
    """
    in_allowed_ontologies = {'MFO': False , 'BPO':False, 'CCO':False}
    if use_has_go:
        if not has_go(sp_rec):
            # print "not has go",sp_rec.entry_name
            return in_allowed_ontologies,[]
    go_list = parse_go(sp_rec)
    for go_rec in go_list:
        if go_rec['evidence'] in allowed_ec:
            in_allowed_ontologies[go_rec['ontology']] = True
    return in_allowed_ontologies, go_list
    

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

def species_go_filter(sp_handle, taxon_id,allowed_ec,go_ec_filter,midfix="_"):
    """
    sp_handle: handle to swissprot file
    taxon_id: NCBI taxon ID
    allowed_ec: list of allowed evidence codes
    go_ec_filter: function, either exclusive or inclusive
    """
    target_id = int(taxon_id+"0000001")
    outhandle = {}
    outhandle['MFO'] = open("sp_species.%s.MFO.%s.tfa" %
                            (taxon_id,midfix),"w")
    outhandle['BPO'] = open("sp_species.%s.BPO.%s.tfa" % 
                            (taxon_id,midfix),"w")
    outhandle['CCO'] = open("sp_species.%s.CCO.%s.tfa" % 
                            (taxon_id,midfix),"w")
    outhandle_all = open("sp_species.%s.all.%s.tfa" % 
                         (taxon_id,midfix) ,"w")
    for inrec in sp.parse(sp_handle):
        if taxon_id in inrec.taxonomy_id:
            in_allowed = go_ec_filter(inrec,allowed_ec)[0]
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

def species_go_data(sp_handle, taxon_id, allowed_ec,
                   go_ec_filter,midfix="_",use_has_go=True):
    """
    Write only GO data for a given species. Use the GO filter function
    provided in the arguments do determine which swissprot entries will be 
    written.
    """
    outhandle = {}
    outhandle['MFO'] = open("go_species.%s.MFO.%s.go" %
                            (taxon_id,midfix),"w")
    outhandle['BPO'] = open("go_species.%s.BPO.%s.go" % 
                            (taxon_id,midfix),"w")
    outhandle['CCO'] = open("go_species.%s.CCO.%s.go" % 
                            (taxon_id,midfix),"w")
    for sp_rec in sp.parse(sp_handle):
        if taxon_id not in sp_rec.taxonomy_id:
            continue
        if use_has_go and not has_go(sp_rec):
            continue
        in_allowed, go_list = go_ec_filter(sp_rec,allowed_ec)
        # print sp_rec.entry_name, in_allowed
        if use_has_go and not go_list:
            continue
        for onto in in_allowed:
            if in_allowed[onto]:
                for go_rec in go_list:
                    if go_rec['evidence'] not in allowed_ec:
                        continue
                    if go_rec['ontology'] != onto:
                        continue

                    outhandle[onto].write("%s\t%s\t%s\n" % 
                        (sp_rec.entry_name,go_rec['go_id'], go_rec['evidence']))
    for i in outhandle:
        outhandle[i].close()

def sp2cafa(sp_go_handle, sp2cafa_handle,out_handle):
    sp2cafa_dict = {}
    for inrec in sp2cafa_handle:
        sp_id, cafa_id = inrec.strip().split()
        sp2cafa_dict[sp_id] = cafa_id
    for inrec in sp_go_handle:
        sp_id, go_id = inrec.strip().split()[:2]
        if sp_id not in sp2cafa_dict:
#            print sp_id
            continue
        else:
            cafa_id = sp2cafa_dict[sp_id]
        
        out_handle.write("%s\t%s\n" % (cafa_id, go_id))

def all_sp2cafa(sp2cafa_path):
    for inpath in glob.glob("go_species.*.go"):
        g1, tax_id, onto, midfix, g2 = inpath.split(".")
        outpath = ("%s.%s.%s.%s.cafa" % (g1, tax_id, onto, midfix))
        sp2cafa(open(inpath), open(sp2cafa_path), open(outpath,"w"))

        

def make_sp_sql_table(sp_handle, out_handle, taxa_ids=[]):
    """
    Make a table from swissprot that includes the following fields:
    SP_ID, GO, Ontology, evidence code
    """
    for sp_rec in sp.parse(sp_handle):
        go_list = parse_go(sp_rec)
        if not go_list:
            continue
        if taxa_ids and (not (sp_rec.taxonomy_id[0] in taxa_ids)):
            continue
        for go_entry in go_list:
            out_handle.write("%s\t%s\t%s\t%s\n" % (sp_rec.entry_name,
                           go_entry['go_id'], 
                           go_entry['ontology'],
                           go_entry['evidence']))

                          
                        
def call_make_sp_sql_table(sp_handle, out_handle, taxa_id_handle):

    taxa_ids = [inrec.split()[0].strip() for inrec in taxa_id_handle]
    print taxa_ids
    make_sp_sql_table(sp_handle, out_handle, taxa_ids)
        

if __name__ == '__main__':

#    species_filter(gzip.open(sys.argv[1]), sys.argv[2])
#    species_go_data(open(sys.argv[1]),sys.argv[2],
#                      exp_ec, go_ec_inclusive_filter,"inclusive-EXP")
#    species_go_data(open(sys.argv[1]),sys.argv[2],
#                      exp_ec, go_ec_inclusive_filter,"inclusive-EXP")
#    all_sp2cafa(sys.argv[1])
    call_make_sp_sql_table(open(sys.argv[1]), open(sys.argv[2],"w"),
                           open(sys.argv[3]))
