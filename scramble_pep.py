#!/usr/bin/env python
'''
Created on 15 Mar 2011

Functions to generate scrambled control peptides.
Peptides are represented as a string of letters, corresponding to the
 20 standard amino acids.
Lower-case letters represent d-Amino acids, and upper-case represents the 
standard l-Amino acids.

@author: Fergal
'''
from __future__ import with_statement #For python 2.5
import itertools
try:
    from itertools import permutations
except ImportError, e:
    def permutations(iterable, r=None):
        # permutations('ABCD', 2) --> AB AC AD BA BC BD CA CB CD DA DB DC
        # permutations(range(3)) --> 012 021 102 120 201 210
        pool = tuple(iterable)
        n = len(pool)
        r = n if r is None else r
        if r > n:
            return
        indices = range(n)
        cycles = range(n, n-r, -1)
        yield tuple(pool[i] for i in indices[:r])
        while n:
            for i in reversed(range(r)):
                cycles[i] -= 1
                if cycles[i] == 0:
                    indices[i:] = indices[i+1:] + indices[i:i+1]
                    cycles[i] = n - i
                else:
                    j = cycles[i]
                    indices[i], indices[-j] = indices[-j], indices[i]
                    yield tuple(pool[i] for i in indices[:r])
                    break
            else:
                return

import re
import sys
import os.path
#import string
import random
try:
    from math import factorial
except ImportError: #python 2.5...
    def factorial(num):
        out = 0
        for i in [n+1 for n in range(num)]:
                if out == 0:
                        out = i
                else:
                        out *= i
        return out

import warnings
#Hide slimfinder deprecation warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import csv
import cStringIO
import traceback
from pprint import pprint
import operator
import optparse

#Non-built in modules
import peptide_occurance
#try:
from suds.client import Client
from PepLibGen.StructGen import synthrules
#except ImportError, e:
#    sys.stderr.write("Could not find required module. ")
#    sys.stderr.write(e.message)

    
#Rich Edwards slimfinder modules
from slimsuite import rje_slim
from slimsuite import rje_slimlist

#################### Module Constants #########################################

#ELMS - use global so that elms aren't read from file multiple times.
elms = []

#Propery swap dictionary, for generating good control peptides. 
property_swap_dict = {'V':'R',
                 'T':'K',
                 'I':'H',
                 'A':'D',
                 'G':'E',
                 'F':'P',
                 'W':'S',
                 'C':'M',
                 'N':'L',
                 'Q':'Y',
                 'R':'V',
                 'K':'T',
                 'H':'I',
                 'D':'A',
                 'E':'G',
                 'P':'F',
                 'S':'W',
                 'M':'C',
                 'L':'N',
                 'Y':'Q'}

#Charged amino acids
charged_residues = ['H', 'R', 'K', 'D', 'E']
#Positively charged amino acids
pos_charged_residues = ['H','R','K']
strict_pos_charged_residues = ['K','R']
#Negatively charged amino acids
neg_charged_residues = ['E','D']
#Uncharged
uncharged_residues = ['I','F','L','W','A','M','P','C','N','V','G','S','Q','Y','T']
#Non-polar residues
nonpolar_residues = ['F','A','L','M','I','W','P','V']
#Small amino acids
small_residues = ['G','A','S','C','P','T','D','N','V']
#hydrophobic residues
hydrophobic_residues = ['C','A','G','I','L','V','T','M','F','Y','W']
#Branched amino-acids
branched_residues = ['V', 'I', 'L']
#AMino-acids with large protecting groups (Fmoc synthesis)
large_protgroup_residues = ['H','R','C','Q','N']

#Filename where short linear motifs will be stored.
#In the module directory...
elm_file = os.path.join(os.path.dirname(__file__), "elm_regex_definitions")

#Cap on number of  randomly generated peptides: to prevent combinatorial
#explosion
MAX_SCRAMBLED_PEPTIDES = 30

MAX_TRIES = 1000

################### Generation Functions #####################################


def _get_pairs(sequence):
    '''
    Returns lists of pairs of peptides that make up a sequence.
    Returns 2 lists, one starting from the start, one offset one residue from
    the start. This means there will be one or two 'unpaired' residues. 

    >>> _get_pairs("FERGAL") 
    ['FE', 'RG', 'AL'], ['F', 'ER', 'GA', 'L']
    >>> _get_pairs("ASDFGHI")
    ['AS', 'DF', 'GH', 'I'], ['A','SD','FG','HI']
    
    '''
    pair1 = [sequence[i:i+2] for i in range(len(sequence)) if not i % 2]
    pair2 = [sequence[0]] +\
     [sequence[i+1:i+3] for i in range(len(sequence)) if not i % 2 and sequence[i+1:i+3]]
    return pair1, pair2

def generate_scrambled_sequences(sequence, pairs=False, filter_funcs=[],
                                 max_out=MAX_SCRAMBLED_PEPTIDES):
    '''
    Generates control peptides by scrambling the sequence.
    
    Can scramble the sequence directly, or use adjacent pairs of peptides
    if pairs=True.
    Accepts functions to filter and modify functions, by setting filter_funcs 
    to a list of functions that accept a peptide sequence as argument and return
    True or False to filter the peptide.
    Accepts a modification function, by setting the mod_func argument,
     which should accept a sequence as an 
    argument, and returna modified sequence.
    '''
    if pairs:
        sequence_pairs = _get_pairs(sequence)        
    #If the sequence is short, directly return a generator for all of them
    #Arbritrary choice of value - factorial(7) is 5040, anything above this 
    #"large enough" that randomly shuffling the sequence will not give 
    #duplicates for a long time.
    #Shorter than that we will generate the whole set 
    #Note: set() will remove duplicates    
    if not pairs and len(sequence) < 8:
        all_controls =  [''.join(control) for control 
                         in set(permutations(sequence, len(sequence))) 
                         if ''.join(control) != sequence]
        for i in range(MAX_SCRAMBLED_PEPTIDES):
            if all_controls:
                yield all_controls.pop(random.randrange(len(all_controls)))
            else:
                break
    else:
        #Shuffle, and return new control peptides
        unsucessful = 0
        already_yielded = []
        while len(already_yielded) < MAX_SCRAMBLED_PEPTIDES:
            if pairs:
                #Use list() to take qcopy of one pair - will be modified in place
                candidate_control = list(random.choice(sequence_pairs))
            else:
                candidate_control = list(sequence)
            #In place list shuffle
            random.shuffle(candidate_control)
            candidate_control = ''.join(candidate_control)
            if candidate_control != sequence: 
                if candidate_control not in already_yielded:
                    already_yielded.append(candidate_control)
                    yield candidate_control
                    unsucessful = 0
            else:
                unsucessful += 1
                if unsucessful > MAX_TRIES:
                    break

def generate_random_chargedist_control(original_seq, opt='keep', filters = [],
                                  max_out=MAX_SCRAMBLED_PEPTIDES):
    '''
    Generates control peptides based on the orginal peptide.
    Takes into account the charge distribution: set the opt parameter to 
    'keep' to return control of the same distribution, or set it to 'opposite' 
    to get the opposite.
    '''
    options = {'keep':(pos_charged_residues, 
                       neg_charged_residues, uncharged_residues),
               'opposite':(neg_charged_residues, 
                           pos_charged_residues, uncharged_residues),
               'uncharged_to_alanine':(pos_charged_residues, 
                                       neg_charged_residues, ['A'])}
    peps1, peps2, peps3 = options[opt]
    returned = []
    unsucessful = 0
    while len(returned) < max_out:
        new_seq = ''
        for resi in original_seq:
            if resi in pos_charged_residues:
                new_seq += random.choice(peps1)
            elif resi in neg_charged_residues:
                new_seq += random.choice(peps2)
            else:
                new_seq += random.choice(peps3)
        #Make sure that the new_seq passes all filters
        if not [f(new_seq) for f in filters if not f(new_seq)]:
            if new_seq != original_seq and new_seq not in returned:
                returned.append(new_seq)
                yield new_seq
                unsucessful = 0
            else:
                unsucessful += 1
                if unsucessful > MAX_TRIES:
                    break
 
def generate_scrambled_chargedist_controls(sequence, opt='keep',
                                           max_out=MAX_SCRAMBLED_PEPTIDES):
    uncharged = [resi for resi in sequence 
                 if resi not in strict_pos_charged_residues 
                 and resi not in neg_charged_residues]
    new_seq = ''
    #substitutions = {'keep':{'neg':'E','pos':'K'},
    #                 'opposite':{'neg':'K','pos':'E'}}
    out_list = []
    tries = 0
    while len(out_list) < max_out:  
        new_seq = ''
        #Copy and shuffle the list
        resi_set = list(uncharged)
        random.shuffle(resi_set)
        for resi in sequence: 
            if resi in neg_charged_residues:
                #new_seq += substitutions[opt]['neg']
                if opt == 'keep':
                    new_seq += resi
                else:
                    new_seq += strict_pos_charged_residues[0]
            elif resi in strict_pos_charged_residues:
                #new_seq += substitutions[opt]['pos']
                if opt == 'keep':
                    new_seq += resi
                else:
                    new_seq += neg_charged_residues[0]
            else:
                new_seq += resi_set.pop()
        if new_seq not in out_list and new_seq != sequence:
            out_list.append(new_seq)
            yield new_seq     
        else:
            tries += 1
            if tries > MAX_TRIES: 
                #Arbitrary max times to try and generate a new output sequence.
                #If an original output cannto be found, exit.
                break

def generate_alanine_scans(sequence):
    for i in range(len(sequence)):
        alanine_scan = sequence[:i] + 'A' + sequence[i+1:]
        if alanine_scan != sequence:
            yield alanine_scan
                
#Modification functions
                      
def property_swap(sequence, strict=True):
    '''
    Takes in a sequence, and swaps amino-acids by differing properties.
    If strict=True, will raise an exception if non-standard amino acids are used 
    in the sequence.
    if strict=False, non standard amino-acids will be left where they are.
    
    >>> property_swap('KVGFFKR')
    'TREPPTV'
    >>> property_swap('FERGALDUFFY')
    Traceback (most recent call last):
        ...
    KeyError: 'U'
    >>> property_swap('FERGALDUFFY', strict=False)
    'PGVEDNAUPPQ'

    '''
    try:
        return ''.join([property_swap_dict[resi] for resi in sequence])
    except KeyError:
        #Probably lower case letters in the sequence, work around.
        new_sequence = ''
        for resi in sequence:
            if resi in property_swap_dict:
                new_sequence += property_swap_dict[resi]
            elif resi.upper() in property_swap_dict:
                new_sequence += property_swap_dict[resi.upper()].lower()
            elif strict:
                #Don't reraise old error, as resi will be different
                raise KeyError(resi)
            else:
                new_sequence += resi
        return new_sequence
      
def reverse(sequence):
    '''
    Reverses the sequence
    '''
    if sequence[::-1] != sequence:
        return sequence[::-1]   

def reversed_pairs(sequence):
    '''
    Returns two sets of reversed pairs of residues.
    '''
    p1, p2 = _get_pairs(sequence)
    return [s for s in [''.join(p1[::-1]), ''.join(p2[::-1])] if s != sequence]

def opposite_charge(sequence):
    '''
    Replaces D, E with K, R; K with E.
    ''' 
    new_seq = ''
    for resi in sequence:
        if resi in ['D','E']:
            new_seq += 'K'
        elif resi in ['R', 'K']:
            new_seq += 'E'
        else:
            new_seq += resi
    return new_seq

        
##################### Filtering functions ####################################


def get_amino_locations(sequence, aminoset):
    '''
    Returns a list of positions of all <aminoset> residues in <sequence>
    '''
    locs = []
    for resi in aminoset:
        location = 0
        while True:
            location = sequence.find(resi, location)
            if location != -1:
                locs.append(location)
                location += 1
            else:
                break
    locs.sort()
    return locs

def well_spaced(sequence, aminoset):
    '''
    Returns True if residues in <aminoset> are well spaced in <sequence>.
    Well distributed means at least a gap of 1 residue between each aminoset 
    residue, as well as aminoset residues being spaced at least 25% of the 
    average distance apart.
    '''
    locations = get_amino_locations(sequence, aminoset)
    if not locations:
        return True
    spacings = [locations[i+1]-locations[i] for i in range(len(locations)-1)]
    expected_spacing = len(sequence)/float(len(locations))
    #Arbitrary - good spacing means all residues are seperated by at least 50%
    #of the average spacing, or over one residue
    good_spacings = [s for s in spacings if s > 0.25*expected_spacing and s > 1]
    if good_spacings == spacings:
        return True
    else:
        return False
        
def biased_to_terminus(sequence, aminoset):
    '''
    Returns 0 for roughly unbiased, -1 for bias to start (nterm) and 1 for bias 
    towards end (cterm)
    '''
    retval = 0
    locations = get_amino_locations(sequence, aminoset)
    #Split sequence in two. If length if not even, must include middle residues
    #in both sides, for balance.
    if not len(sequence) % 2:
        start, end = sequence[:len(sequence)/2], sequence[len(sequence)/2:]
    else:
        start, end = sequence[:(len(sequence)/2)+1], sequence[len(sequence)/2:]
    start_count = sum([start.count(resi) for resi in aminoset])
    end_count = sum([end.count(resi) for resi in aminoset])
    if (start_count > end_count*2 and end_count) or start_count > end_count+2:
        retval = -1
    elif (end_count > start_count*2 and start_count) or end_count > start_count+2:
        retval = 1
    return retval

def set_elm_regexes():
    '''
    Returns the list of short linear motif regexes.
    Will try and get them from a file in the current directory, if that doesn't
    exist, will read them from the ELM database, and attempt to save them
    to a text file for future use, one regex per line.
     
    >>> set_elm_regexes()
    True
     
    '''
    global elms
    read_from_file =False
    if os.path.exists(elm_file):
        #get definitions from file
        with open(elm_file, 'r') as f:
            lines=f.readlines()
            if lines:
                read_from_file = True
                #ELMs file format:
                #1 ELM per line
                # Name, regex, IC
                elms = []
                for line in lines:
                    line = line.split()
                    elms.append((line[0].strip().strip(','), 
                                line[1].strip().strip(','), 
                                float(line[2].strip().strip(','))))
    if not read_from_file:
        # get definitions from server
        elms = get_elm_regexes_from_server()
        #print elms
        #Write to file
        try:
            with open(elm_file, 'w') as f:
                for triple in elms:
                    triple = (triple[0],
                              triple[1],
                              str(triple[2]))
                    f.write(', '.join(triple)+"\n")
        except Exception, e:
            print "Could not write ELMs to file", e
            os.remove(elm_file)
    if elms:
        return True
    else:
        return False

def get_elm_regexes_from_server():
    '''
    Returns all short linear motif regexes from the ELM database.
    '''
    wsdl = 'http://elm.eu.org/webservice/ELMdb.wsdl'
    client = Client(wsdl)
    #Return list of regexes.
    all_elms = client.service.getAllELMs()
    return [(elm['Identifier'], 
             elm['Regex'], 
             get_motif_information_content(elm['Regex'])) for elm in all_elms]

def difficult_peptide(sequence):
    '''
    Returns a string describing the difficult synthesis features of a peptide.
    
    '''
    difficult = synthrules.test_motifs(sequence)
    if not difficult:
        difficult = []
    if not synthrules.test_charge(sequence):
        difficult.append("Does not have a charged residue every 5 residues")
    if not well_spaced(sequence, branched_residues):
        difficult.append("Branched residues are not well spaced")
    if biased_to_terminus(sequence, branched_residues) == 1:
        difficult.append("Several branched amino-acids close to C-terminus")
    if not well_spaced(sequence, large_protgroup_residues):
        difficult.append("Contains residues with large Fmoc protecting groups close together")
    return difficult

def get_matching_motifs(sequence):
    '''
    Returns a list of motif names that match the sequence
    '''
    motifs = []
    if not elms:
        print "Retrieving ELM definitions."
        set_elm_regexes()
        #pprint(elms)
    for name, motif, _ in elms:
        if re.search(motif, sequence):
            motifs.append(name)
    return motifs

def get_motif_information_content(regex):
    '''
    Calculates the information content of a regular expression/SLiM.
    Uses slimsearch modules to convert the slim to PRESTO format (made up of 
    individual SLiM elements, seperated by '-', and evaluates and sums the
    information content of each SLiM element to return a value.
    For regular expression containing '|' that map to several possible PRESTO
    SliMs, returns the lowest found information content.
    '''
    #Cleanup
    #slimFromPattern does not handle '?' in regex properly so:
    regex = regex.replace('?', '{0,1}')
    #Doesn't like '$' either, but they don't change the information content so:
    #Bit of a corner case...
    regex = regex.replace('|$', '')
    regex = regex.replace('$', '') 
    if '|' in regex and not try_split(regex):
        #Need to split regex into all possibilites, as PRESTO format cannot
        #handle converting the '|' in the regex.
        #SLiMlist class can handle this, the _splitMotif method will break up 
        #the regex and instantiate rje_slim.SLiM() objects for each regex, 
        #which can be accessed from the slimlist.list['Motif'] variable.
        #Also, need to wrap the regex in '()' for splitting, as the 
        #_splitMotif method expects '|' to be inside '()' at some stage
        regex = '(' + regex + ')'
        slimlist = rje_slimlist.SLiMList()
        slimlist._splitMotif('',regex,'')
        presto_slims = [slim.info['Slim'] for slim in slimlist.list['Motif']]
    elif '|' in regex:
        #slimlist seems to have trouble splitting regexes that have no '('
        #groups, which can be easily split by try_split
        presto_slims =\
         [rje_slim.slimFromPattern(regex) for regex in try_split(regex)]
    else:
        #The .slimFromPattern function derives the PRESTO slim representation
        #directly from the regex.
        presto_slims = [rje_slim.slimFromPattern(regex)]
    #pprint(presto_slims)
    #Function calculates the information content for a presto-format SLiM.
    ic = lambda slim: sum([rje_slim.elementIC(ele) for ele in slim.split('-')])
    ics = [ic(slim) for slim in presto_slims]
    #pprint(ics)
    return min([ic(slim) for slim in presto_slims])
    
def try_split(regex):
    '''
    Simplistic algorithm for splitting a regular expression by '|'.
    Only works if there is exactly 1 '|' and the '|' is not inside a group
    '''
    if not regex.count('|') == 1:
        return False
    regexes = regex.split('|')
    for regex in regexes:
        #Check that both regexes are properly parenthesised
        depth = 0
        for char in regex:
            if char == '(':
                depth += 1
            elif char == ')':
                depth -= 1
        if depth != 0:
            return False
    return regexes
    
def reassemble_groups(groups):
    return ''.join(['(' + group + ')' for group in groups])

############################ Main and associated functions ####################

def yield_controls(sequence, synth_filter=False, 
                   motif_filter=False, chatty=False):
    '''
    Main function: outputs a .csv with columns of control peptides under a 
    title: e.g. scrambled, pair scrambled, etc.
    ''' 
    # 26/05/2011: TODO List
    # Single peptide controls: reverse, opposite charge, reversed pairs
    # Same / Opposite charge controls should use original amino acids.
    # Pos/neg controls should use only Ks and Es (No H, D, R)
    # 
    if chatty:
        printout = sys.stdout
    else:
        printout = open(os.devnull, 'w')
    #Main dictionary for controls: will be a dictionary with the name of the 
    #control type as a key, and a generator of values as a control.
    #def generate_scrambled_sequences(sequence, pairs=False, filter_funcs=[], 
    #                                   mod_func='', 
    #                                   max=MAX_SCRAMBLED_PEPTIDES):
    base_filters = []
#        if synth_filter:
#            base_filters.append(synthesisable)
#        if motif_filter:
#            base_filters.append(no_motifs)
    all_controls = {}
    all_controls['Original_sequence'] = [sequence]
    
    print>>printout, "Generating scrambled controls."
    all_controls['Randomly_scrambled'] = \
    list(generate_scrambled_sequences(sequence, filter_funcs=base_filters))
    
    all_controls['Reversed_sequence']\
     = [reverse(sequence)] if reverse(sequence) else []
    # all_controls['Opposite_charge_distribution'] = \
    # [opposite_charge(sequence)]
     
    print>>printout, "Generating pair scrambled controls."
    all_controls['Pairs_scrambled'] =\
    list(generate_scrambled_sequences(sequence, pairs=True, 
                                      filter_funcs=base_filters))
    
    all_controls['Reversed_pairs'] = reversed_pairs(sequence)
    
    print>>printout, "Generating same-charge scrambled peptides"
    all_controls["Fixed_charge_and_scrambled"] = \
    list(generate_scrambled_chargedist_controls(sequence))
    
    print>>printout, "Generating Alanine-scan scrambled peptides"
    all_controls["Alanine_scans"] = list(generate_alanine_scans(sequence))
     
    print>>printout, "Generating random same-charge controls."
    all_controls["Random_peptide_with_same_charge"] = \
    list(generate_random_chargedist_control(sequence, filters=base_filters))
            
    #CSV output
    #Write headings, but not if filtering has happened(do that later)
    #Write headings for csv file..
    headings = ['Control Type:',
                'Sequence:', 
                'logP(o/w):',
                'Synthesis problems (lower is better):', 
                'Synthesis summary:', 
                'Matching ELM motifs',
                'Matching ELM motif with highest informativeness',
                'Highest ELM Informativeness',
                'Largest peptide shared with original sequence',
                'Length of longest peptide shared with original sequence',
                'Needleman-wunsch alignment score to original sequence',
                'Levenshtien distance to original sequence']
    yield headings
    #Write the original sequence
    #writer.writerow(['Original_sequence',all_controls['Original_sequence']])
    #pprint(all_controls)
    #pprint(all_controls.items())
    sort_order = ['Original_sequence', 
                  'Reversed_sequence', 
                  #'Opposite_charge_distribution', 
                  'Reversed_pairs',
                  'Randomly_scrambled', 
                  'Pairs_scrambled', 
                  'Fixed_charge_and_scrambled',
                  'Alanine_scans',
                  #'Scrambled_opposite_charge',
                  'Random_peptide_with_same_charge', 
                 # 'Random_opposite_charge_distribution']
                  ]
    sorted_controls = sorted(all_controls.items(), 
                             key = lambda controls: sort_order.index(controls[0]))
    if not elms:
        set_elm_regexes()
    for key, values in sorted_controls:
        for i, peptide in enumerate(values):
            matching_motifs = get_matching_motifs(peptide)
            sorted_matching_motifs = sorted([elm for elm in elms if elm[0] in matching_motifs], 
                                        key=operator.itemgetter(2))
            highest_elm_ic_name = sorted_matching_motifs[0][0] if sorted_matching_motifs else ''
            highest_elm_ic_value = sorted_matching_motifs[0][2] if sorted_matching_motifs else ''
            row = [key + '_' + str(i+1) if len(values) > 1 else key, 
                   peptide, 
                   round(synthrules.test_hydrophobicity(peptide), 2),
                   len(difficult_peptide(peptide)) if difficult_peptide(peptide) else 0,
                   ', '.join(difficult_peptide(peptide)),
                   ', '.join(matching_motifs),
                   highest_elm_ic_name,
                   highest_elm_ic_value,
                   peptide_occurance.largest_shared_peptide(peptide, sequence),
                   len(peptide_occurance.largest_shared_peptide(peptide, sequence)),
                   peptide_occurance.needleman_wunsch(peptide, sequence)[0],
                   peptide_occurance.levenshtein_distance(peptide, sequence)]
            yield row     


def write_controls(sequence, synth_filter=False, 
                   motif_filter=False, outfile = True, outname=False, chatty=True):
    '''
    Main function: outputs a .csv with columns of control peptides under a 
    title: e.g. scrambled, pair scrambled, etc.
    ''' 
    # 26/05/2011: TODO List
    # Single peptide controls: reverse, opposite charge, reversed pairs
    # Same / Opposite charge controls should use original amino acids.
    # Pos/neg controls should use only Ks and Es (No H, D, R)
    # 
    #Output will be .csv format
    if not outfile:
        out = cStringIO.StringIO()
    else:
        name = outname if outname else sequence+'_controls.csv'
    if chatty:
        printout = sys.stdout
    else:
        printout = open(os.devnull, 'w')
    try:
        if outfile:
            try:
                out = open(name, 'wb')
            except IOError: #File is open in another applicaiton
                out = open(sequence + '_controls_1.csv', 'wb')
        
        writer = csv.writer(out)
        for row in yield_controls(sequence, synth_filter=synth_filter, 
                                  motif_filter=motif_filter, chatty=chatty):
            writer.writerow(row)
        writer.writerow(row)     
        if not outfile:
            return out.getvalue()
        else:
            return True
    finally:
        if not printout is sys.stdout:
            printout.close()
        #Cleanup open files
        try:
            out.close()
        except:
            print "Error cleaning up..."
            pass
        
def main():
    usage = "Usage: scramble_pep.py [options] input_sequence(s)"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option('-o', '--outfile-name', type='string', default="",
                      help="""Set name of output file. 
                              Default is <input_sequence>_controls.csv. 
                              Option can only be set for one input sequence""")
    options, args = parser.parse_args()
    if len(args) != 1 and options.outfile_name != "":
        print "-o --outfile-name option only valid for a single input sequence"
    elif len(args) == 1:
        print "Generating controls for", args[0]
        outname = options.outfile_name if options.outfile_name else True
        write_controls(args[0], outname=outname)
    else:
        for arg in args:
            print "Generating controls for", arg
            write_controls(arg)           
            print
    print "Done."

############# Where the magic happens... #####################################

if __name__ == '__main__': 
    main()
    
    
