#!/usr/bin/env python
import sys
import glob
from collections import defaultdict
from Bio.UniProt import GOA as upg
GO_EXP_EC = set(['EXP','IDA','IPI','IMP','IGI','IEP'])

def filter_in_IEA(handle):
    outhandle = open(handle.name+".IEA","w")
    outhandle.write('!gaf-version: 2.0\n')
    for inrec in upg.gafiterator(handle):
        if inrec['Evidence'] == 'IEA':
            upg.writerec(inrec,outhandle) 
    outhandle.close()

def filter_in_experimental(handle):
    outhandle = open(handle.name+".exp_evidence","w")
    outhandle.write('!gaf-version: 2.0\n')
    for inrec in upg.gafiterator(handle):
        if upg.record_has(inrec,{'Evidence':  GO_EXP_EC}):
            upg.writerec(inrec,outhandle) 
    outhandle.close()


def became_experimental(handle_iea, handle_exp):
    """Identify Electronically annotated proteins that became experimentally
    annotated later
    handle_iea: gaf file with proteins exclusively IEA annotated
    handle_exp: gaf file with proteins experimentally annotated
    Note: files should contain only one ontology (either MFO, BPO or
    CCO) to make sense.
    """
    expdict = {}
    outhandle = open(handle_exp.name+".became_exp","w")
    outhandle.write('!gaf-version: 2.0\n')
    # First read experimental into memory
    for exprec in upg.gafbyproteiniterator(handle_exp):
        expdict[exprec[0]['DB_Object_ID']] = exprec
    # Now read in the non-experimental
    for iearec in upg.gafbyproteiniterator(handle_iea):
        prot_id = iearec[0]['DB_Object_ID']
        if prot_id in expdict:
            upg.writebyproteinrec(expdict[prot_id],outhandle)
    outhandle.close()
     
        
def extract_taxa(handle, taxalist):
    """
    Create a GAF file from multiple taxa
    taxalist is a list of strings of taxid. Don't use list of int
    """
    outfiles = {}
    header = "!gaf-version: 2.0\n"
    for taxid in taxalist:
        outfiles[taxid] = open("%s.taxon.%s" % (handle.name,taxid),'w')
        outfiles[taxid].write(header)
    for inrec in upg.gafiterator(handle):
        cur_taxid = inrec['Taxon_ID'][0].split(':')[1]
        if cur_taxid in taxalist:
            upg.writerec(inrec, outfiles[cur_taxid])
    for i in outfiles:
        outfiles[i].close()
        

def extract_taxon(handle,in_taxid):
    """
    Create a GAF file from a single taxon
    """
    header = "!gaf-version: 2.0\n"
    if isinstance(in_taxid, int):
        taxid = str(in_taxid)
    taxid = in_taxid.strip()
    outfile = open("%s.taxon.%s" % (handle.name,taxid),'w')
    outfile.write(header)
    for inrec in upg.gafiterator(handle):
        if inrec['Taxon_ID'][0].split(':')[1] == taxid:
            upg.writerec(inrec, outfile)
    outfile.close()

def split_to_ontologies(handle):
    """Splits a GAF file into three ontology files
    """
    header = "!gaf-version: 2.0\n"
    out_mfo = open("%s.MFO" % handle.name,'w')
    out_bpo = open("%s.BPO" % handle.name,'w')
    out_cco = open("%s.CCO" % handle.name,'w')
    out_bpo.write(header)
    out_mfo.write(header)
    out_cco.write(header)
    for inrec in upg.gafiterator(handle):
        if inrec['Aspect'] == 'F':
            upg.writerec(inrec,out_mfo)
        elif inrec['Aspect'] == 'P':
            upg.writerec(inrec,out_bpo)
        elif inrec['Aspect'] == 'C':
            upg.writerec(inrec,out_cco)
        else:
            raise ValueError, 'unknown ontology aspect %s' % inrec['Aspect']
    out_mfo.close()
    out_bpo.close()
    out_cco.close()

def species_stats(handle):
    """Statistics for species distributions in a gaf file"""
    taxa_count = defaultdict(int)
    for prot_rec in upg.gafbyproteiniterator(handle):
        taxa_count[prot_rec[0]['Taxon_ID'][0]] += 1
    return taxa_count

def print_species_stats(taxa_count,outhandle):
    taxa_list = [(i[1],i[0]) for i in list(taxa_count.items())]
    taxa_list.sort()
    taxa_list.reverse()
    taxa_names = read_taxa_names(open('names.dmp'))
    for count,taxon_id in taxa_list:
        try: 
            taxon_name = taxa_names[taxon_id.split(':')[1]]
        except KeyError:
            taxon_name = taxon_id
        if count > 1:
            outhandle.write("%s\t%d\n" % (taxon_name, count))
    outhandle.close()
def get_species_stats(handle):
    taxa_count = species_stats(handle)
    print_species_stats(taxa_count,open(handle.name+".taxa_stats","w"))

def read_taxa_names(names_handle):
    taxa_names = {}
    for inrec in names_handle:
        taxon_id, taxon_name = inrec.split('|')[:2]
        taxa_names[taxon_id.strip()] = taxon_name.strip()
    return taxa_names
    
def exclusive_IEA(goa_reclist):
    exclusive = True
    for rec in goa_reclist:
        if rec['Evidence'] != 'IEA':
            exclusive = False
            break
        continue 
    return exclusive

def has_experimental(goa_reclist):
    retval = False
    for rec in goa_reclist:
        if upg.record_has(rec,{'Evidence': GO_EXP_EC}):
            retval = True
            break
    return retval

def all_hasnt_experimental(handle):
    outhandle = open(handle.name+".noexp","w")
    outhandle.write('!gaf-version: 2.0\n')

    for protrec in upg.gafbyproteiniterator(handle):
        if not has_experimental(protrec):
            for outrec in protrec:
                upg.writerec(outrec, outhandle)
    outhandle.close()

def all_exclusive_IEA(handle):
    outhandle = open(handle.name+".exclusive_IEA","w")
    outhandle.write('!gaf-version: 2.0\n')

    for protrec in upg.gafbyproteiniterator(handle):
        if exclusive_IEA(protrec):
            for outrec in protrec:
                upg.writerec(outrec, outhandle)
    outhandle.close()


if __name__ == '__main__':
    # typical target prep
    inlist = [i.split()[0].strip() for i in  open(sys.argv[2]).readlines()]
    extract_taxa(open(sys.argv[1]),inlist)
    sys.exit()
#    extract_taxon(open(sys.argv[1]), sys.argv[2])
#    split_to_ontologies(open(sys.argv[1]))
#    for infile in glob.iglob("*.[108.MFO|108.BPO|108.CCO]"):
    for infile in ['gene_association.goa_uniprot.BPO',
                   'gene_association.goa_uniprot.MFO',
                   'gene_association.goa_uniprot.CCO']:
        print infile
#        all_hasnt_experimental(open(infile))
        filter_in_experimental(open(infile))
#        became_experimental(open(sys.argv[1]+".108.exclusive_IEA"), open(sys.argv[1]))
#    get_species_stats(open(sys.argv[1]))
#    taxa_count = species_stats(open(sys.argv[1]))
#    print_species_stats(taxa_count,open(sys.argv[1]+".taxa_count","w"))
    pass    

