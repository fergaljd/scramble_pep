#!/usr/bin/env python
'''
Parse and pickle/unpickle the UNIPROT reference proteome.
'''
from __future__ import with_statement
import re
import os
import cPickle

human_proteins_fasta = 'uniprot-organism-human-reference-proteome.fasta'
human_proteins_pickle = 'uniprot-organism-human-reference-proteome.pkl'

def protein_supplier(fasta):
    '''
    yields protein sequences and annotations, given a UNIPROT fasta file.
    '''
    with open(fasta) as proteins:
        protein = {}
        annotations = ['uniprot_id', 'name', 'organism', 'GN', 'PE', 'SV']
        for line in proteins:
            if line[0] == '>':
                if protein: 
                    yield protein
                protein = {}
                header = line[1:].strip()
                #Header line looks like
                #>sp|Q6P3R8|NEK5_HUMAN Serine/threonine-protein kinase Nek5 OS=Homo sapiens GN=NEK5 PE=1 SV=1
                protein['db'], protein['accession'], refs = header.split('|')
                ref_pattern = re.compile(r"""(?P<uniprot_id>\w+)\s
                                             (?P<name>.+)\s
                                             (?P<organism>OS=(.+)\s)? 
                                             (?P<GN>GN=(.+)\s)? 
                                             (?P<PE>PE=(.+)\s)? 
                                             (?P<SV>SV=(.+))?""", re.VERBOSE)
                match = ref_pattern.match(refs)
                for annotation in annotations:
                    protein[annotation] = match.group(annotation) 
            else:
                if protein.has_key('sequence'):
                    protein['sequence'] += line.strip()
                else:
                    protein['sequence'] = line.strip()
        yield protein
        
def pickle_proteome(fasta=human_proteins_fasta, 
                    pickle_file=human_proteins_pickle):
    '''
    Parses and pickles the proteome FASTA file.
    '''
    proteins = list(protein_supplier(fasta))
    with open(pickle_file, 'wb') as picklefile:
        cPickle.dump(proteins, picklefile, protocol=2)
    return proteins
        
def load_proteome(pickle_file=human_proteins_pickle):
    '''
    Loads the pickled human proteome.
    '''
    with open(pickle_file, 'rb') as pickle:
        return cPickle.load(pickle)
    
def get_proteins(fasta = human_proteins_fasta, 
                 picklefile=human_proteins_pickle):
    if os.path.exists(picklefile):
        return load_proteome(picklefile)
    else:
        return pickle_proteome(fasta)
