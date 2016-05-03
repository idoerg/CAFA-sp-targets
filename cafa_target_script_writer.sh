
source query_species.sh | gawk 'BEGIN {FS=","} {print "python filter_sp_species.py /home/idoerg/work/CAFA3/uniprot_sprot.dat " $2 " &\t# " $1}' >| create_cafa_targets.sh
