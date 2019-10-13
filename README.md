CAFA-sp-targets
===============

Generates targets for CAFA from a SwissProt file. 

Prerequesites: biopython http://biopython.org

Download the latest SwissProt file from here: 
(ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz)

Targets for CAFA 4 will be provided from the September 18, 2019 release of UniprotKB/SwissProt (2019_08)
(ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2019_08/knowledgebase/knowledgebase2019_08.tar.gz)

This downloads a large (15GB) tarball containing a TrEMBL file and a SwissProt file. Only the Swissprot file will be used.

To run:
```
$ python cafa_sp_species.py swissprot_file taxid
```
Where tax_id is the NCBI taxonomy id of the taxon for which you wish to generate CAFA targets. E.g. 9606 is Homo sapiens.

See the NCBI Taxonomy ID page for more details: http://www.ncbi.nlm.nih.gov/taxonomy

Output:

1) sp_species.[taxon_id].tfa

2) sp_species.[taxon_id].MFO.noexp.tfa

3) sp_species.[taxon_id].BPO.noexp.tfa

4) sp_species.[taxon_id].CCO.noexp.tfa 

5) sp_species.[taxon_id].all.noexp.tfa


The files contain:
1) FASTA files with all the proteins in [taxon_id] (replaced by an actual number, e.g. sp_species.9606.tfa for human.)

2-4) FASTA files for proteins that are not experimentally annotated in the MFO, BPO or CCO ontologies.

5) FASTA file containing proteins that are not experimentally annotated in *any* GO ontology.


