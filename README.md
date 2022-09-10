# fetch_afdb_pymol

after sourcing the script on pymol 

``run afdb.py`` or ``run afdb_noBS.py``

then type and run

``af UNIPROTID``

This fetches the corresponding alphafold-database structure model.

For example,

``af A0A452S449``

# Dependency

afdb.py requires  BeautifulSoup to show the Uniprot fasta information.

afdb_noBS.py is a reduced version which does not require BautifulSoup but does not show Uniprot information (just fetches model)


This script is modification of the ``fetch`` command from PyMOL source code.

There're many unused pieces of code, so it's welcomed to cleaning them up on PR.
