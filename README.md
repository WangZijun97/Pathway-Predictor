Pathway-Predictor

Simple python script to predict biological pathways based on protein interaction data from Reactome's database.

Use: Download, run `python "src_cached/main.py" start_acc end_acc iterations`

`start_acc` and `end_acc` refer to Accession numbers of proteins as used by Uniprot
Max distance to explore = 2 x `iterations`

Please download Protein Atlas Database (proteinatlas.org) and STRING human database (string-db.org), and store them within `/data`.
Ensure that the Protein Atlas Database is `proteinatlas.json` and the STRING human database is `9606.protein.links.v11.0.txt`.

`/data/conversions.json` is cached data from uniprot.org's mapping service between identifiers.
