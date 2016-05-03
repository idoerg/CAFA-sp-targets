#!/usr/bin/env python

import os
import sys
from collections import defaultdict
import sqlite3
#import GOAParser # In future, this module will be replaced with the Bio-Python module GOA.py in Bio.UniProt
from Bio.UniProt import GOA as GOAParser # In future, this module will be replaced with the Bio-Python module GOA.py in Bio.UniProt

'''
           This script extracts all swiss prot proteins , for a particular taxon, 
           along with its associated information from uniprot-goa files.It takes as 
           input a gene-association file, gene protein information file and a taxon 
           id as input. The first 2 files ideally should be of the same version to 
           avoid missing any proteins and they can be downloaded from uniprot-goa ftp site 
           ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/. 

           Usage :
           python create_SP_only_protein_dataset.py <gene association file> 
                                                    <gene protein information file> <taxon id>

           Output:
           .gaf file contains all information in gaf format and 
           .db file that puts all information into a sqlite3 database table


'''

def parse_gpi(infile, taxon=''):
    
    sp_id = defaultdict()

    infile_handle = open(infile, 'r')
    parser = GOAParser.gpi_iterator(infile_handle)

    for rec in parser:
        print rec.keys()
        if not rec.has_key('Gene_Product_Properties'):
            print "This version of the gp information file does not contain all required information"
            sys.exit(1)
        else:
            break

    for rec in parser:
        taxid = rec['Taxon'].split(':')[1].strip()
        db = rec['Gene_Product_Properties'][0].split('=')[1].strip()
        if db.startswith('Swiss-Prot') and taxon == taxid:
            sp_id[rec['DB_Object_ID']] = 1

    return sp_id



def insert_into_db(record, taxon, GAFFIELDS):

    conn = None
    try:
        conn = sqlite3.connect("filtered_gaf_file_with_only_sp_ids_for_" + taxon + ".db")    
        c = conn.cursor()
        c.execute("drop table if exists SP_gaf")

        if len(GAFFIELDS) == 17:
            c.execute("create table SP_gaf(db varchar(20), db_id varchar(20), db_symbol varchar(20), " + \
                          "qualifier varchar(20), GO_Term varchar(40), db_ref varchar(40), " + \
                          "evidence varchar(3), With varchar(20), ontology char(1), db_name varchar(100), " + \
                          "synonym varchar(200), db_type varchar(40), taxid varchar(20), Date date, " + \
                          "source varchar(40), ann_ext varchar(20), gene_product varchar(20))")

            c.executemany("insert into SP_gaf values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", record)
            conn.commit()
        elif len(GAFFIELDS) == 15:
            c.execute("create table SP_gaf(db varchar(20), db_id varchar(20), db_symbol varchar(20), " + \
                          "qualifier varchar(20), GO_Term varchar(40), db_ref varchar(40), " + \
                          "evidence varchar(3), With varchar(20), ontology char(1), db_name varchar(100), " + \
                          "synonym varchar(200), db_type varchar(40), taxid varchar(20), Date date, source varchar(40))")
            c.executemany("insert into SP_gaf values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", record)
            conn.commit()

    except sqlite3.Error as e:
        print "Error %s:" % e.args[0]
        sys.exit(1)

    finally:
        if conn:
            conn.close()



def extract_gaf(rec, outfile, GAFFIELDS, record, sp_id, taxon):

    t = ()

    if sp_id.has_key(rec['DB_Object_ID']):    
        GOAParser.writerec(rec,outfile,GAFFIELDS)

        if len(GAFFIELDS) == 15:
            t = (rec['DB'], rec['DB_Object_ID'], rec['DB_Object_Symbol'], ('|'.join(rec['Qualifier'])), 
                 rec['GO_ID'], ('|'.join(rec['DB:Reference'])), rec['Evidence'], ('|'.join(rec['With'])), 
                 rec['Aspect'], rec['DB_Object_Name'], ('|'.join(rec['Synonym'])), rec['DB_Object_Type'], 
                 ('|'.join(rec['Taxon_ID'])), rec['Date'], rec['Assigned_By'])

        elif len(GAFFIELDS) == 17:
            t = (rec['DB'], rec['DB_Object_ID'], rec['DB_Object_Symbol'], ('|'.join(rec['Qualifier'])), 
                 rec['GO_ID'], ('|'.join(rec['DB:Reference'])), rec['Evidence'], ('|'.join(rec['With'])), 
                 rec['Aspect'], ('|'.join(rec['DB_Object_Name'])), ('|'.join(rec['Synonym'])), 
                 rec['DB_Object_Type'], ('|'.join(rec['Taxon_ID'])), rec['Date'], rec['Assigned_By'], 
                 rec['Annotation_Extension'], rec['Gene_Product_Form_ID'])
        
        record.append(t)

    return record


if __name__ == '__main__':

    gaf_file = sys.argv[1]
    gpi_file = sys.argv[2]
    taxon = sys.argv[3]

    outfile = open("filtered_gaf_file_with_only_sp_ids_for_" + taxon + ".gaf" ,'w')
    gaf_handle = open(gaf_file, 'r')
    record = []

    sp_id = parse_gpi(gpi_file, taxon)
    parser = GOAParser.gafiterator(gaf_handle)

    for rec in parser:
        if len(rec) == 15:
            GAFFIELDS = GOAParser.GAF10FIELDS
            break
        elif len(rec) == 17:
            GAFFIELDS = GOAParser.GAF20FIELDS
            break

    for rec in parser:
        record = extract_gaf(rec, outfile, GAFFIELDS, record, sp_id, taxon)
    
    new_record = tuple(record)
    insert_into_db(new_record, taxon, GAFFIELDS)
