CAFA-sp-targets
===============

CAFA-sp-targets.py

Generates targets for CAFA from a SwissProt file. 

Prerequesites: biopython http://biopython.org
Download the latest SwissProt file from here: 
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

Targets for CAFA 2 will be provided from the June 2013 releae of SwissProt:
ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2013_06/knowledgebase/knowledgebase2013_06.tar.gz

This downloads a large (15GB) tarball containing a TrEMBL file and a SwissProt file. Only the Swissprot file will be used.

To run:

$ python cafa_sp_targets.py swissprot_file taxid

Where swissprot_file
