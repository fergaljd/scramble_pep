# scramble pep


Command line script for generating scrambled peptides: scramble_pep.py
Can generate random peptides based on an original sequence, or peptides with opposite/similar charge and hydrophilicity properties.
Filters out peptides with known binding motifs (based on short linear motifs: SLiMs/ELMs)
Can also filter out synthetically difficult peptides.

Requires installation of PepLibGen, available in my CycloPs reposory. 
Also requires the installation of the SUDS library in order to query short linear motid definitions.
