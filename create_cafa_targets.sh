#!/bin/bash
export CAFASOFT='/home/idoerg/soft/CAFA-sp-targets-py3'
export CAFADATA='/home/idoerg/work/CAFA4'
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  9606 &	# "homo sapiens[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  10090 &	# "mus musculus[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  10116 &	# "Rattus norvegicus"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  3702 &	# "Arabidopsis thaliana[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  83333 &	# "Escherichia coli K-12[all names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  7227 &	# "Drosophila melanogaster[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  287 &	# "Pseudomonas aeruginosa[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  559292 &	# "Saccharomyces cerevisiae ATCC 204508"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  559292 &	# "Schizosaccharomyces pombe ATCC 24843"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  7955 &	# "Danio rerio[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  44689 &	# "Dictyostelium discoideum[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  243273 &	# "Mycoplasma genitalium ATCC 33530"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  6239 &	# "Caenorhabditis elegans[All Names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  226900 &	# "Bacillus cereus ATCC 14579"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  4577 &	# "Zea Mays [All names]"
python $CAFASOFT/filter_sp_species.py $CAFADATA/uniprot_sprot.dat  9823 &	# "Sus scrofa" 
