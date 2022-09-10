# fetch_afdb_pymol

after sourcing the script for pymol (run afdb.py)

``af UNIPROTID``

fetches the corresponding alphafold-database structure model.

For example,

``af A0A452S449``

# dependdency

afdb.py requires  BeautifulSoup to show the Uniprot fasta information.

afdb_noBS.py is reduced version which does not require BautifulSoup but does not show Uniprot information (just fetches model)


This script is modification of the ``fetch`` command from PyMOL source code.

There're many unused pieces of code, so it's welcomed to cleaning them up on PR.
