#!/usr/bin/env python
'''
Calculated Peptide Occurance statistics for the reference human proteome
'''
import sys
import re
import pdb
import traceback
from timeit import Timer
import cProfile

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo
import proteome_pickler
       
def get_subpeptides(sequence):
    '''
    Breaks a peptide into all possible sub-peptides
    >>> get_subpeptides('ASD')
    >>> ['ASD', 'AS', 'SD', 'A', 'S', 'D']
    '''
    subpeptides = []
    window_size = len(sequence)
    while window_size > 0:
        if window_size == len(sequence):
            subpeptides.append(sequence)
            window_size -= 1
        else:
            offset = 0
            while (window_size + offset) <= len(sequence):
                subpeptides.append(sequence[offset:window_size+offset])
                offset += 1
            window_size -= 1
    return subpeptides

def calculate_peptide_occurance(sequence):
    '''
    Find how often a peptide and its subpeptides are seen in the human proteome. 
    '''
    subpeptides = get_subpeptides(sequence)
    occurances = [0]*len(sequence)
    for protein in proteome_pickler.get_proteins():
        matches = [len(pep) for pep in subpeptides 
                   if pep in protein['sequence']]
        for m in matches:
            occurances[m] += 1
    return occurances

def largest_shared_peptide(s1, s2):
    '''
    Find the largest subpeptide shared between 2 peptides.
    '''
    try:
        return max([subpep for subpep in get_subpeptides(s1) 
                    if subpep in get_subpeptides(s2)], key=len)
    except ValueError: #No shared subpeptides,max of empty sequence
        return ""
        
def levenshtein_distance(s1, s2):
    '''
    Calculates the minimum distance between 2 strings, with insertion, deletion
    and substitution being allowed. 
    Better than Haming distance, which only allows subsitution.
    '''
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if not s1:
        return len(s2)
    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def needleman_wunsch(s1, s2):
    '''
    Uses BioPythons local alignment code to implememt the needleman-wunsch
    algorithm. Uses the Blosum62 scoring matrix, and no gap penalty.
    Returns the score, and the s1 and s2 alignments.
    '''
    #Smith-waterman sets negative scoring 
    #Set gap creation penalty to 0.
    aln_s1, aln_s2, score, begin, end =\
     pairwise2.align.globaldx(s1, s2, MatrixInfo.blosum62)[0]
    return score, aln_s1, aln_s2
     
            
if __name__ == '__main__':
    try:
        #t = Timer("calculate_peptide_occurance(sys.argv[1])", 
        #          "from __main__ import calculate_peptide_occurance")
        #print t.timeit(number=50)
        cProfile.run('calculate_peptide_occurance(sys.argv[1])')
    except:
        type, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)

