# fetch_afdb_pymol

After sourcing the script on pymol by ``run afdb.py`` or ``run afdb_noBS.py``

Then type and run

``af UNIPROTID``

This fetches the corresponding alphafold-database structure model.

For example,

``af A0A452S449``

fetches Free fatty acid receptor 2

Similarly, after sourcing the script on pymol by ``run esmatlas.py```

Then type and run

``esm MGnifyID``

This fetches the corresponding ESM-atlas structure model.

For example,

``esm MGYP002537940442``

will work. 

# Dependency

afdb.py requires  BeautifulSoup to show the Uniprot fasta information.

afdb_noBS.py is a reduced version which does not require BautifulSoup but does not show Uniprot information (just fetches model)

# Notes

This script is modification of the ``fetch`` command from PyMOL source code.

There're many unused pieces of code inside, so it's welcomed to cleaning them up on PR.
