#!/usr/bin/python

# rje_sequence - DNA/Protein sequence object
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_sequence
Description:  DNA/Protein sequence object
Version:      1.16
Last Edit:    16/11/10
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the Sequence Object used to store sequence data for all PEAT applications that used DNA or
    protein sequences. It has no standalone functionality.

    This modules contains all the methods for parsing out sequence information, including species and source database,
    based on the format of the input sequences. If using a consistent but custom format for fasta description lines,
    please contact me and I can add it to the list of formats currently recognised.

Uses general modules: copy, os, random, re, sre_constants, string, sys, time
Uses RJE modules: rje, rje_disorder
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, os, random, re, sre_constants, string, sys, time
#########################################################################################################################
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_disorder
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Separated Sequence object from rje_seq.py
    # 1.1 - Rudimentary opt['GeneSpAcc'] added
    # 1.2 - Modified RegExp for sequence detail extraction
    # 1.3 - Added list of secondary accession numbers and hasID() method to check ID and all AccNum (and gnspacc combos)
    # 1.4 - Added Peptide Design methods
    # 1.5 - Added storing of case in a dictionary self.dict['Case'] = {'Upper':[(start,stop)],'Lower':[(start,stop)]}
    # 1.6 - Added disorder and case masking
    # 1.7 - Added FudgeFactor and AA codes
    # 1.8 - Added position-specific AA masking
    # 1.9 - Added EST translation functions. Fixed fudging. Added dna() method.
    # 1.10- Fixed sequence name bug
    # 1.11- Added recognition of UniRef
    # 1.12- Added AA masking
    # 1.13- Added Taxonomy list and UniProt dictionary for UniProt sourced sequences (primarily).
    # 1.14- Added maskRegion()
    # 1.15- Added disorder proportion calculations.
    # 1.16- Added additional Genbank and EnsEMBL BioMart sequence header recognition.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add descriptions of recognised formats to documentation.
    '''
#########################################################################################################################
### Regular Expressions for Sequence Details Extraction
re_acconly = re.compile('^(\S+)\s*$')
re_plain = re.compile('^(\S+)\s+(\S.*)$')
re_gnspec = re.compile('^(\S+)_(\S+)')
re_gnspacc = re.compile('^(\S+)_([A-Za-z0-9]+)\/(\S+)')
re_gn_sp__acc = re.compile('^(\S+)_([A-Za-z0-9]+)__(\S+)')
re_uniprot = re.compile('^(sp|tr|sw|uniprot)\|(\S+)\|(\S+)_(\S+)')
re_uniref = re.compile('^(UniRef\S+)_(\S+)\s.*\Tax=(\S.+)\sRepID=(\S+)')
re_uniprot2 = re.compile('^(\S+)_(\S+)\s+\((\S+)\)')
re_uniprot3 = re.compile('^([A-Za-z0-9\-]+)\|([A-Za-z0-9]+)_([A-Za-z0-9]+)\s+')
re_disprot = re.compile('^DisProt\|(DP\d+)\|')
re_genbank = re.compile('^(\S+)\|(\d+)\|')
re_ncbi = re.compile('^gi\|(\d+)\|(\S+)\|(\S+)\|')
re_ncbi2 = re.compile('^gi\|(\d+)\|(\S+)\|(\S+)')
re_unigene = re.compile('^\S+\|UG\|(\S+)')
re_ensemblpep = re.compile('^(\S+)\s+pep:(\S+)')
re_ipi = re.compile('^>*IPI:(IPI\d+\.\d+)\|.+Tax_Id=(\d+)')
re_hprd = re.compile('^ID_HPRD_(\d+)')
re_hprd2 = re.compile('^(\d+)\|(\d+_\d+)\|(\S+)\|')
re_pipe = re.compile('^(\S+)_(\S+)\|(\S+)')
re_db_pipe = re.compile('^(\S+)\|(\S+)')
re_tigr = re.compile('^(\S+)\s.+\{(\S.+)\}')
re_flybase = re.compile('^(\S+) type=.+name=(\S+);.+species=(\S+);')
re_jgi = re.compile('^jgi\|(\S+)\|(\d+)\|(\S+)')
re_enst = re.compile('^(ENS\S+)\|ENS(\S+)\|')
ipi_taxa = {'3702':'ARATH', '7955':'BRARE', '9031':'CHICK', '9606':'HUMAN', '10090':'MOUSE', '10116':'RAT'}
#########################################################################################################################
genetic_code = {'UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'S','UCA':'S','UCG':'S','UAU':'Y','UAC':'Y','UAA':'*','UAG':'*','UGU':'C','UGC':'C','UGA':'*','UGG':'W',
                'CUU':'L','CUC':'L','CUA':'L','CUG':'L','CCU':'P','CCC':'P','CCA':'P','CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','CGU':'R','CGC':'R','CGA':'R','CGG':'R',
                'AUU':'I','AUC':'I','AUA':'I','AUG':'M','ACU':'T','ACC':'T','ACA':'T','ACG':'T','AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R','AGG':'R',
                'GUU':'V','GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A','GCG':'A','GAU':'D','GAC':'D','GAA':'E','GAG':'E','GGU':'G','GGC':'G','GGA':'G','GGG':'G'}
aa_code_3 = {'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His','I':'Ile','L':'Leu','M':'Met',
             'P':'Pro','Q':'Gln','R':'Arg','S':'Ser','T':'Thr','V':'Val','X':'Unk','Y':'Tyr','*':'STOP'}
aa_3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'MET': 'M', 'THR': 'T', 'PRO': 'P', 'STOP': '*', 'HIS': 'H',
           'PHE': 'F', 'ALA': 'A', 'GLY': 'G', 'ILE': 'I', 'LEU': 'L', 'ARG': 'R', 'UNK': 'X', 'VAL': 'V', 'GLU': 'E',
           'TYR': 'Y'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Sequence Class for Individual sequences                                                                 #
#########################################################################################################################
class Sequence(rje.RJE_Object):
    '''
    Individual Sequence Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of sequence
    - Type = Type of sequences
    - Sequence = Actual sequence
    - Description = Description of sequence list (if desired)
    - ID = Sequence identifier (unique 'word')
    - AccNum = Accession Number
    - DBase = Source database (for AccNum)
    - Gene = Gene Symbol
    - Species = species name
    - SpecCode = SwissProt Species Code
    - Format = Name Format
    - NCBI = NCBI GenBank gi number
    
    Stat:numeric

    Opt:boolean
    - RevComp = Whether sequence has been reverse complemented or not [False]

    List:list
    - IsDisordered = List of True/False values for each residue, whether disordered or not
    - Secondary ID = List of secondary IDs and Accession numbers
    - Taxonomy = List of taxonomic levels for source species
    
    Dict:dictionary
    - Case = Stores case of original input sequence as tuples {'Upper':[(start,stop)],'Lower':[(start,stop)]}
    - UniDAT = UniProt data dictionary (if sequence read from UniProt entry)

    Obj:RJE_Objects
    - Disorder = rje_disorder.Disorder object containing disorder prediction results. (Must run self.disorder() first!)
    '''
    ### Attributes
    def seqLen(self): return len(self.info['Sequence'])
    def aaLen(self): return len(self.info['Sequence']) - string.count(self.info['Sequence'],'-')
    def aaNum(self): return self.aaLen()
    def nonX(self):
        if self.dna(): return self.nonN()
        return self.aaLen() - string.count(self.info['Sequence'].upper(),'X')
    def nonN(self): return self.aaLen() - string.count(self.info['Sequence'].upper(),'N')
    def dna(self): return self.seqType() in ['RNA','DNA']   # Whether a DNA (or RNA) sequence
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### <a> ### Basics 
        self.infolist = ['Name','Type','Description','Sequence','ID','AccNum','DBase','Gene','Species','SpecCode','Format',
                         'PreMask','MaskSeq','NCBI']
        self.statlist = []
        self.optlist = ['RevComp']
        self.listlist = ['IsDisordered','Secondary ID','Taxonomy']
        self.dictlist = ['Case','UniDAT']
        self.objlist = ['Disorder']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'SpecCode':'UNK','NCBI':''})
        self.dict['Case'] = {'Upper':[],'Lower':[]}
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)
                self._cmdReadList(cmd,'opt',['Yeast'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        return
#########################################################################################################################
    ### <2> ### Sequence Data (from name) extraction and updates                                                        #
#########################################################################################################################
    def addSequence(self,sequence='',case=True,caselist=[],stripnum=False): ### Adds sequence to object, extracting case information
        '''
        Adds sequence to object, extracting case information.
        >> sequence:str [''] = sequence to add to self.info['Sequence']
        >> case:bool [True] = whether to store case information in self.dict['Case']
        >> caselist:list [] = list of case positions to over-ride read in case (self.list['Case'] from rje_seq)
        >> stripnum:boolean [False] = whether to strip numbers from sequences
        '''
        try:
            ### Strip Space and Numbers ###
            sequence = string.join(string.split(sequence),'')
            if stripnum:
                for badness in ['0','1','2','3','4','5','6','7','8','9']: sequence = string.replace(sequence,badness,'')
            ### Update Info ###
            self.info['Sequence'] = sequence.upper()
            ### Case Dictionary ###
            self.dict['Case'] = {'Upper':[],'Lower':[]}
            if caselist:
                switch = {'Upper':'Lower','Lower':'Upper'}
                seq = {'Upper':sequence.upper(),'Lower':sequence.lower()}
                clist = [0]
                for i in caselist:
                    try: clist.append(string.atoi(i))
                    except: self.log.errorLog('Cannot use "%s" from caselist!' % i)
                c = 'Upper'
                while clist:
                    r = clist.pop(0)
                    sequence = sequence[:r] + seq[c][r:]
                    c = switch[c]
            if case: self.dict['Case'] = caseDict(sequence)
        except: self.log.errorLog('Problem with rje_sequence.addSequence()')
#########################################################################################################################
    def getSequence(self,case=False,gaps=True):   ### Returns sequence
        '''
        Returns sequence.
        >> case:bool [False] = whether to use self.dict['Case'] information.
        >> gaps:bool [True] = whether to leave gaps in the sequence
        '''
        sequence = self.info['Sequence'].upper()
        if case:
            for (start,stop) in self.dict['Case']['Lower']:
                sequence = sequence[:start] + sequence[start:(stop+1)].lower() + sequence[(stop+1):]
        if gaps: return sequence
        else: return string.replace(self.info['Sequence'],'-','')
#########################################################################################################################
    def reverseComplement(self,rna=False):  ### Converts to reverse complement of sequence (upper case)
        '''Converts to reverse complement of sequence (upper case).'''
        self.info['Sequence'] = reverseComplement(self.info['Sequence'],rna=rna)
        self.opt['RevComp'] = not self.opt['RevComp']
#########################################################################################################################
    def trimPolyA(self):    ### Removes 3' As
        '''Removes 3' As.'''
        while self.info['Sequence'][-1:] == 'A': self.info['Sequence'] = self.info['Sequence'][:-1]
#########################################################################################################################
    def extractDetails(self,gnspacc=False):   ### Extracts ID, Acc# etc. from sequence name
        '''
        Extracts ID, Acc# etc. from sequence name.
        '''
        ### <a> ### General Case
        try:
            self.info['FullName'] = self.info['Name']
            self.info['Name'] = string.replace(self.info['Name'],'||','|')  #!# Dodgy!
            if re_plain.match(self.info['Name']):
                match = re_plain.match(self.info['Name'])
                self.info['Format'] = "plain"
                grp = match.groups()
                self.info['ID'] = grp[0]
                self.info['AccNum'] = grp[0]
                self.info['Description'] = grp[1]
                if rje.matchExp('gi:(\d+)',self.info['Description']): self.info['NCBI'] = rje.matchExp('gi:(\d+)',self.info['Description'])[0]
            elif re_acconly.match(self.info['Name']):
                match = re_acconly.match(self.info['Name'])
                self.info['Format'] = "acconly"
                grp = match.groups()
                self.info['ID'] = grp[0]
                self.info['AccNum'] = grp[0]
                #print self.info['Name'],self.info['AccNum']
                #raw_input()
            else:
                self.log.errorLog('Problem finding sequence ID etc. from %s' % self.info['Name'])
                return 0
        except:
            self.log.errorLog('Error extracting basic sequence name information.')
        ### <b> ### gene_SPECIES/Acc# and other Formats
        try:
            if re_gnspacc.match(self.info['Name']):
                self.info['Format'] = 'gnspacc'
                match = re_gnspacc.match(self.info['Name'])
                grp = match.groups()
                self.info['Gene'] = grp[0]
                self.info['SpecCode'] = grp[1]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = grp[2]
                if self.info['Gene'].upper() == self.info['Gene']:
                    if self.info['Gene'] != self.info['AccNum']: self.info['DBase'] = 'sprot'
                    else: self.info['DBase'] = 'trembl'
            elif re_enst.match(self.info['Name']):
                details = string.split(self.info['Name'],'|')   # re_jgi = re.compile('^jgi\|Thaps3\|(\d+)\|(\S+)')
                self.info['DBase'] = self.info['Format'] = 'ENST'
                self.info['Gene'] = details[-1].lower()
                self.info['AccNum'] = details[0]
                self.info['Description'] = details[2]
                self.info['EnsG'] = details[1]
                self.info['SpecCode'] = details[0][3:6]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                #self.deBug(self.info['ID'])
            elif re_jgi.match(self.info['Name']):
                details = string.split(self.info['Name'],'|')   # re_jgi = re.compile('^jgi\|Thaps3\|(\d+)\|(\S+)')
                self.info['Gene'] = 'jgi'
                if details[1] == 'Thaps3':
                    self.info['SpecCode'] = 'THAPS'
                    self.info['Species'] = 'Thalassiosira pseudonana'
                    if details[1][-2:] == 'bd': self.info['AccNum'] = 'TPSBD%s' % details[2]
                    else: self.info['AccNum'] = 'TPSJGI%s' % details[2]
                elif details[1] == 'Emihu1':
                    self.info['SpecCode'] = 'EMIHU'
                    self.info['Species'] = 'Emiliania huxleyi'
                    self.info['AccNum'] = 'EHUXJGI%s' % details[2]
                elif details[1] == 'Nemve1':
                    self.info['SpecCode'] = 'NEMVE'
                    self.info['Species'] = 'Nematostella vectensis'
                    self.info['AccNum'] = 'NEMVEJGI%s' % details[2]
                else:
                    self.info['SpecCode'] = details[1].upper()
                    self.info['AccNum'] = 'JGI%s' % details[2]
                self.info['Description'] = details[3]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                #self.deBug(self.info['ID'])
            elif re_gnspacc.match(self.info['Description']):
                self.info['Format'] = 'gnspacc fragment'
                match = re_gnspacc.match(self.info['Description'])
                grp = match.groups()
                self.info['Gene'] = grp[0]
                self.info['SpecCode'] = grp[1]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = grp[2]
                if self.info['Gene'].upper() == self.info['Gene']:
                    if self.info['Gene'] != self.info['AccNum']: self.info['DBase'] = 'sprot'
                    else: self.info['DBase'] = 'trembl'
            elif re_disprot.match(self.info['Name']):
                self.info['Format'] = 'DisProt'
                self.info['DBase'] = 'disprot'
                self.info['AccNum'] = re_disprot.match(self.info['Name'])[0]
                self.info['Gene'] = 'dp'
            elif re_gn_sp__acc.match(self.info['Name']):
                self.info['Format'] = 'gn_sp__acc'
                match = re_gn_sp__acc.match(self.info['Name'])
                grp = match.groups()
                self.info['Gene'] = grp[0]
                self.info['SpecCode'] = grp[1]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = grp[2]
                if self.info['Gene'].upper() == self.info['Gene']:
                    if self.info['Gene'] != self.info['AccNum']: self.info['DBase'] = 'sprot'
                    else: self.info['DBase'] = 'trembl'
            elif re_hprd.match(self.info['Name']):
                self.info['DBase'] = self.info['Format'] = 'HPRD'
                self.info['ID'] = self.info['Name']
                self.info['AccNum'] = 'HPRD%s' % rje.regExp(re_hprd,self.info['Name'])[0]
                self.info['Gene'] = 'hprd'
                self.info['SpecCode'] = 'HUMAN'
            elif re_hprd2.match(self.info['Name']):
                self.setInfo({'DBase':'HPRD','Format':'HPRD_Raw','Gene':'hprd','SpecCode':'HUMAN'})
                (self.info['ID'],sv,self.info['AccNum']) = rje.regExp(re_hprd2,self.info['Name'])
            elif re_uniprot.match(self.info['Name']):
                self.info['Format'] = 'uniprot'
                detail = rje.regExp(re_uniprot,self.info['Name'])
                self.info['Gene'] = detail[2]
                self.info['SpecCode'] = detail[3]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = detail[1]
                if self.info['Gene'] == self.info['AccNum']:
                    self.info['DBase'] = 'trembl'
                    self.info['Gene'] = 'tr'
                else:
                    self.info['DBase'] = 'sprot'
            elif re_uniref.match(self.info['Name']):    # re.compile('^(UniRef\S+)_(\S+)\s.*\Tax=(\S.+)\sRepID=(\S+)')
                self.info['Format'] = 'uniref'
                detail = rje.regExp(re_uniref,self.info['Name'])
                if detail[3].find('_') > 0: [self.info['Gene'],self.info['SpecCode']] = string.split(detail[3],'_')
                else:
                    self.info['Gene'] = 'uref90'
                    self.info['SpecCode'] = getSpecCode(detail[2])
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = detail[1]
                self.info['DBase'] = detail[0]
            elif re_uniprot2.match(self.info['Name']):
                self.info['Format'] = 'uniprot2'
                detail = rje.regExp(re_uniprot2,self.info['Name'])
                self.info['Gene'] = detail[0]
                self.info['SpecCode'] = detail[1]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = detail[2]
                if self.info['Gene'] == self.info['AccNum']:
                    self.info['DBase'] = 'trembl'
                    self.info['Gene'] = 'tr'
                else:
                    self.info['DBase'] = 'sprot'
            elif re_uniprot3.match(self.info['Name']):
                self.info['Format'] = 'uniprot3'
                detail = rje.regExp(re_uniprot3,self.info['Name'])
                self.info['Gene'] = detail[1]
                self.info['SpecCode'] = detail[2]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = detail[0]
                if self.info['Gene'] == self.info['AccNum']:
                    self.info['DBase'] = 'trembl'
                    self.info['Gene'] = 'tr'
                else:
                    self.info['DBase'] = 'sprot'
            elif re_ipi.match(self.info['Name']):
                self.info['DBase'] = self.info['Format'] = 'IPI'
                self.info['ID'] = rje.regExp(re_ipi,self.info['Name'])[0]
                self.info['AccNum'] = self.info['ID']
                _taxa = rje.regExp(re_ipi,self.info['Name'])[1]
                if _taxa in ipi_taxa.keys():
                    self.info['SpecCode'] = ipi_taxa[_taxa]
                else:
                    self.info['SpecCode'] = 'UNK'
            elif re_ncbi.match(self.info['Name']):
                self.info['Format'] = 'NCBI'
                [self.info['NCBI'],self.info['DBase'],self.info['AccNum']] = rje.regExp(re_ncbi,self.info['Name'])
                self.info['Gene'] = self.info['DBase']
                if rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',self.info['ID']):
                    [self.info['AccNum'],self.info['Gene'],self.info['SpecCode']] = rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',self.info['ID'])
                    self.info['DBase'] = 'sprot'
                else:
                    self.info['ID'] = '%s_%s' % (self.info['DBase'],self.info['AccNum'])
                    if rje.matchExp('\[(.+)\]\s*$',self.info['Name']):
                        self.info['Species'] = rje.matchExp('\[(.+)\]\s*$',self.info['Name'])[0]
                        while re.search('\[(.+)$',self.info['Species']):
                            self.info['Species'] = rje.matchExp('\[(.+)$',self.info['Species'])[0]
                        self.info['SpecCode'] = getSpecCode(self.info['Species'])
                self.info['Description'] = 'gi:%s %s' % (self.info['NCBI'],self.info['Description'])
            elif re_ncbi2.match(self.info['Name']):
                self.info['Format'] = 'NCBI'
                [self.info['NCBI'],self.info['DBase'],self.info['AccNum']] = rje.regExp(re_ncbi2,self.info['Name'])
                self.info['Gene'] = self.info['DBase']
                self.info['ID'] = '%s_%s' % (self.info['DBase'],self.info['AccNum'])
                if rje.matchExp('\[(.+)\]\s*$',self.info['Name']):
                    self.info['Species'] = rje.matchExp('\[(.+)\]\s*$',self.info['Name'])[0]
                    while re.search('\[(.+)$',self.info['Species']):
                        self.info['Species'] = rje.matchExp('\[(.+)$',self.info['Species'])[0]
                    self.info['SpecCode'] = getSpecCode(self.info['Species'])
                self.info['Description'] = 'gi:%s %s' % (self.info['NCBI'],self.info['Description'])
            elif re_genbank.match(self.info['Name']):
                self.info['Format'] = 'GenBank'
                [self.info['DBase'],self.info['ID']] = rje.regExp(re_genbank,self.info['Name'])
                self.info['AccNum'] = self.info['ID']
            elif re_unigene.match(self.info['Name']):
                self.info['DBase'] = self.info['Format'] = 'UniGene'
                self.info['ID'] = rje.regExp(re_unigene,self.info['Name'])[0]
                self.info['AccNum'] = self.info['ID']
            elif re_flybase.match(self.info['Name']):
                (acc,gene,spec) = rje.regExp(re_flybase,self.info['Name'])
                self.info['DBase'] = 'FlyBase'
                self.info['ID'] = self.info['AccNum'] = acc
                if gene[-3:] in ['-RA','-PA']: gene = gene[:-3]
                self.info['Gene'] = string.split(gene,'\\')[-1]
                self.info['SpecCode'] = '%sRO%s' % (spec[0].upper(),spec[1:3].upper())
            elif re_ensemblpep.match(self.info['Name']):
                self.info['Format'] = 'ensemblpep'
                detail = rje.regExp(re_ensemblpep,self.info['Name'])
                #print detail
                self.info['AccNum'] = detail[0]
                self.info['Gene'] = 'p'
                if detail[1].lower() == 'known': self.info['Gene'] = 'ens'
                elif detail[1].lower() == 'novel': self.info['Gene'] = 'nvl'
                elif detail[1].lower() == 'genscan': self.info['Gene'] = 'scan'
                self.info['SpecCode'] = 'UNK'
                #print self.info['AccNum'], self.info['Gene']
                if self.info['AccNum'].find('ENSAPMP') == 0: self.info['SpecCode'] = 'APIME'
                elif self.info['AccNum'].find('ENSBTAP') == 0: self.info['SpecCode'] = 'BOVIN'
                elif self.info['AccNum'].find('ENSGALP') == 0: self.info['SpecCode'] = 'CHICK'
                elif self.info['AccNum'].find('ENSPTRP') == 0: self.info['SpecCode'] = 'PANTR'
                elif self.info['AccNum'].find('ENSCINP') == 0: self.info['SpecCode'] = 'CIOIN'
                elif self.info['AccNum'].find('ENSCAFP') == 0: self.info['SpecCode'] = 'CANFA'
                elif self.info['AccNum'].find('CG') == 0: self.info['SpecCode'] = 'DROME'
                elif self.info['AccNum'].find('FBpp') == 0: self.info['SpecCode'] = 'DROME'
                elif self.info['AccNum'].find('SINFRUP') == 0: self.info['SpecCode'] = 'FUGRU'
                elif self.info['AccNum'].find('ENSDARP') == 0: self.info['SpecCode'] = 'DANRE'
                elif self.info['AccNum'].find('ENSP0') == 0: self.info['SpecCode'] = 'HUMAN'
                elif self.info['AccNum'].find('ENSANGP') == 0: self.info['SpecCode'] = 'ANOGA'
                elif self.info['AccNum'].find('AGAP') == 0: self.info['SpecCode'] = 'ANOGA'
                elif self.info['Name'].find('group:AMEL') > 0: self.info['SpecCode'] = 'YEAST'
                elif self.info['AccNum'].find('ENSMUSP') == 0: self.info['SpecCode'] = 'MOUSE'
                elif self.info['AccNum'].find('ENSRNOP') == 0: self.info['SpecCode'] = 'RAT'
                elif self.info['AccNum'].find('GSTENP') == 0: self.info['SpecCode'] = 'TETNG'
                elif self.info['AccNum'].find('ENSXETP') == 0: self.info['SpecCode'] = 'XENLA'
                elif self.info['Name'].find('chromosome:SGD') > 0: self.info['SpecCode'] = 'YEAST'
                elif self.info['Name'].find('chromosome:CEL') > 0: self.info['SpecCode'] = 'CAEEL'
                elif self.info['Name'].find('chromosome:WS') > 0: self.info['SpecCode'] = 'CAEEL'
                elif self.info['Name'].find('ENSMMUP') == 0: self.info['SpecCode'] = 'MACMU'
                elif self.info['Name'].find('scaffold:MMUL') > 0: self.info['SpecCode'] = 'MACMU'
                elif self.info['Name'].find('ENSMODP') == 0: self.info['SpecCode'] = 'MONDO'
                elif self.info['Name'].find('scaffold:BROAD') > 0: self.info['SpecCode'] = 'MONDO'
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            elif re_pipe.match(self.info['Name']):
                self.info['Format'] = 'pipe'
                detail = rje.regExp(re_pipe,self.info['Name'])
                self.info['Gene'] = detail[0]
                self.info['SpecCode'] = detail[1]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = detail[2]
                if self.info['Gene'] == self.info['AccNum']:
                    self.info['DBase'] = 'trembl'
                    self.info['Gene'] = 'tr'
                else:
                    self.info['DBase'] = 'sprot'
            elif re_db_pipe.match(self.info['Name']):
                self.info['Format'] = 'db_pipe'
                (self.info['DBase'],self.info['AccNum']) = rje.regExp(re_db_pipe,self.info['Name'])
                self.info['Gene'] = self.info['DBase'].lower()
                self.info['ID'] = self.info['AccNum']
            elif re_gnspec.match(self.info['Name']):
                self.info['Format'] = 'gnspec'
                detail = rje.regExp(re_gnspec,self.info['Name'])
                self.info['Gene'] = detail[0]
                self.info['SpecCode'] = detail[1]
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['AccNum'] = self.info['ID']
                self.info['DBase'] = 'custom'
            elif re_tigr.match(self.info['Name']):
                self.info['DBase'] = self.info['Format'] = 'TIGR'
                (self.info['AccNum'],self.info['Species']) = rje.matchExp(re_tigr,self.info['Name'])
                while rje.matchExp('\S+\|(\S+)',self.info['AccNum']):
                    self.info['AccNum'] = rje.matchExp('\S+\|(\S+)',self.info['AccNum'])[0]
                self.info['Gene'] = 'p'
                self.info['SpecCode'] = getSpecCode(self.info['Species'])
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                
        ### <c> ### EnsEMBL
            if self.info['Gene'] == 'ens':
                self.info['DBase'] = 'ens_known'
            elif self.info['Gene'] == 'nvl':
                self.info['DBase'] = 'ens_novel'
            elif self.info['Gene'] == 'scan':
                self.info['DBase'] = 'ens_scan'
            elif self.info['Gene'] == 'ipi':
                self.info['DBase'] = 'IPI'
            elif self.info['Gene'] == 'ref':
                self.info['DBase'] = 'ens_known'
            #!# Add more Pattern matching and the like for different formats!

            ## Special Kate Database etc. ##
            self.specialDetails()

        ### <d> ### gnspacc
            for g in self.info['Gene'][0:]:
                if not rje.matchExp('([A-Za-z0-9_-])',g) and g not in ['.','#']: self.info['Gene'] = string.replace(self.info['Gene'],g,'')
            if gnspacc:
                if self.info['Gene'].lower() == 'none':
                    self.info['ID'] = '%s_%s' % (self.info['AccNum'],self.info['SpecCode'])
                else:
                    self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
                self.info['Name'] = '%s__%s' % (self.info['ID'],self.info['AccNum'])
                if self.info['Description'].lower() != 'none':
                    self.info['Name'] += ' %s' % self.info['Description']
                self.info['Format'] == 'gn_sp__acc'
            return 1
        except:
            self.log.errorLog('Error extracting advanced sequence name information (%s).' % self.info['Name'])
#########################################################################################################################
    def specialDetails(self):   ### Access special information from Names
        '''Access special information from Names.'''
        try:### ~ [1] ~ Special Kate Database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'].find('TP_') >= 0 and self.info['Name'].find('{Treponema pallidum Nichols}') >= 0:
                self.info['Gene'] == 'p'
                self.info['SpecCode'] = 'TREPA'
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            ### ~ [2] ~ Yeast ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Yeast' in self.opt and self.opt['Yeast']:
                if self.shortName()[4:5] == '_': self.info['AccNum'] = self.shortName()[:4].upper() + self.shortName()[5:]
                if self.shortName()[:4] == 'Kpol': self.setInfo({'SpecCode':'VANPO','Gene':'Kpol'})
                elif self.shortName()[:4] == 'Scas': self.setInfo({'SpecCode':'SACCA','Gene':'Scas'})
                elif self.shortName()[:4] == 'CAGL': self.setInfo({'SpecCode':'CANGA','Gene':'Cgla','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Cgla': self.setInfo({'SpecCode':'CANGA','Gene':'Cgla'})
                elif self.shortName()[:4] == 'Sbay': self.setInfo({'SpecCode':'SACBA','Gene':'Sbay'})
                elif self.shortName()[:4] == 'Anc_': self.setInfo({'SpecCode':'YANC','Gene':'Anc','AccNum':'Anc%s' % self.shortName()[4:]})
                elif self.shortName()[:4] == 'Scer': self.setInfo({'SpecCode':'YEAST','Gene':'Scer'})
                elif rje.matchExp('(Y\S+\d+[WC])',self.shortName()): self.setInfo({'SpecCode':'YEAST','Gene':'Scer','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Agos': self.setInfo({'SpecCode':'ASHGO','Gene':'Agos'})
                elif rje.matchExp('(A\S+\d+[WC])',self.shortName()): self.setInfo({'SpecCode':'ASHGO','Gene':'Agos','AccNum':self.shortName()})
                elif self.shortName()[:3] == 'ZYR': self.setInfo({'SpecCode':'ZYGRO','Gene':'Zrou','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Zrou': self.setInfo({'SpecCode':'ZYGRO','Gene':'Zrou'})
                elif self.shortName()[:4] == 'Klac': self.setInfo({'SpecCode':'KLULA','Gene':'Klac'})
                elif self.shortName()[:4] == 'KLLA': self.setInfo({'SpecCode':'KLULA','Gene':'Klac','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Kwal': self.setInfo({'SpecCode':'KLUWA','Gene':'Kwal'})
                elif self.shortName()[:4] == 'Sklu': self.setInfo({'SpecCode':'SACKL','Gene':'Sklu'})
                elif self.shortName()[:4] == 'SAKL': self.setInfo({'SpecCode':'SACKL','Gene':'Sklu','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'KLTH': self.setInfo({'SpecCode':'LACTH','Gene':'Kthe','AccNum':self.shortName()})
                elif self.shortName()[:4] == 'Kthe': self.setInfo({'SpecCode':'LACTH','Gene':'Kthe'})
                self.info['AccNum'] = string.replace(self.info['AccNum'],'YGOB_','')
                self.info['AccNum'] = string.replace(self.info['AccNum'],'Anc_','Anc-')
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            ### ~ [3] ~ Special E hux EST translation consensi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if (self.shortName()[:3] == 'EHC' or self.shortName()[:4] == 'EHUX') and self.shortName() == self.info['Name']:
                self.info['Gene'] == 'est'
                self.info['SpecCode'] = 'EMIHU'
                self.info['Species'] = 'Emiliania huxleyi'
                self.info['ID'] = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            ### ~ [4] ~ eORF analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            eorf_re = '^(\S+.+)\s(\S+\|.+)gene:(\S+)'
            if rje.matchExp(eorf_re,self.info['Name']):
                eorf = rje.matchExp(eorf_re,self.info['Name'])
                acc = string.split(eorf[0])[0]
                if len(string.split(acc,'_')) > 2: return
                gene = string.split(eorf[1],'|')[0]
                if '_CAEEL' in gene: (gene,spcode) = string.split(gene,'_')
                elif 'WS200' in eorf[1]: spcode = 'CAEEL'
                elif acc[:7] == 'ENSGALT': spcode = 'CHICK'
                elif acc[:7] == 'ENSBTAT': spcode = 'BOVIN'
                elif acc[:4] == 'FBtr': spcode = 'DROME'
                elif acc[:4] == 'ENST': spcode = 'HUMAN'
                elif acc[:7] == 'ENSMUST': spcode = 'MOUSE'
                elif acc[:7] == 'ENSSSCT': spcode = 'PIG'
                elif acc[:7] == 'ENSXETT': spcode = 'XENTR'
                elif acc[:7] == 'ENSDART': spcode = 'DANRE'
                elif 'SGD' in eorf[1]: spcode = 'YEAST'
                if len(gene) < 2 or len(gene) > 6: gene = 'eorf'
                self.info['Gene'] = gene.lower()
                self.info['AccNum'] = acc
                self.info['SpecCode'] = spcode
            elif len(self.info['SpecCode']) == 3:
                if self.info['SpecCode'] == 'T00': (self.info['SpecCode'],self.info['Species']) = ('HUMAN','Homo sapiens')
                
        except: self.errorLog('Special Details problem for %s' % self.info['Name'])
#########################################################################################################################
    def newGene(self,gene='p',keepsp=False,gnspacc=True):     ### Gives sequence new gene
        '''
        Gives sequence new gene.
        >> gene:str = new gene
        >> keepsp:str = whether to keep an UPPER CASE gene identifier
        '''
        if keepsp and (self.info['DBase'] == 'sprot' or self.info['Gene'].find(gene.upper()) >= 0): return  # SPROT            
        self.info['Gene'] = gene
        self.info['ID'] = '%s_%s' % (gene,self.info['SpecCode'])
        if self.info['Format'] == 'gnspacc':
            self.info['Name'] = '%s/%s %s' % (self.info['ID'],self.info['AccNum'],self.info['Description'])
        elif self.info['Format'] in ['gn_sp__acc','uniprot2'] or gnspacc:
            self.info['Name'] = '%s__%s %s' % (self.info['ID'],self.info['AccNum'],self.info['Description'])
        else:
            self.verbose(0,3,'Gene but not Name changed for %s, (%s format).' % (self.shortName(),self.info['Format']),1)
#########################################################################################################################
    def deGap(self):    ### Degaps sequence
        '''Degaps sequence.'''
        self.info['Sequence'] = string.join(string.split(self.info['Sequence'],'-'),'')
#########################################################################################################################
    ### <3> ### Sequence information methods                                                                            #
#########################################################################################################################
    def shortName(self):    ### Returns short name.
        '''Returns short name = first word of name.'''
        try: return string.split(self.info['Name'])[0]
        except:
            self.log.errorLog('Major problem with shortName(%s)' % self.info['Name'])
            raise
#########################################################################################################################
    def seqType(self):  ### Returns (and possible guesses) Sequence Type - Protein/DNA/RNA
        '''
        Returns (and possible guesses) Sequence Type
        - Protein if non-ATGCUN
        - DNA if ATGCN only
        - RNA if AUGCN only
        '''
        if self.info['Type'] == 'None': # Work it out
            if re.search('[DEFHIKLMPQRSVWY]',self.info['Sequence'].upper()): self.info['Type'] = 'Protein'
            elif re.search('U',self.info['Sequence'].upper()): self.info['Type'] = 'RNA'
            else: self.info['Type'] = 'DNA'
        return self.info['Type']            
#########################################################################################################################
    def sameSpec(self,otherseq,unk=False):    ### Returns true if same spec (or SpecCode) or False if not
        '''
        Returns true if same spec (or SpecCode) or False if not.
        >> otherseq:Sequence Object.
        >> unk:bool [False] = value to be returned if either sequence is of unknown species.
        '''
        try:
            for unktest in ['Unknown','UNK']:
                if unktest in [self.info['SpecCode'],otherseq.info['SpecCode']]: return unk
            if self.info['SpecCode'] == otherseq.info['SpecCode']: return True
            return False
        except:
            self.log.errorLog('Major problem with sameSpec()')
            return False
#########################################################################################################################
    def hasID(self,id=None,gnspacc=True,uniprot=True):  ### Checks ID, AccNum and secondary ID for match
        '''
        Checks ID, AccNum and secondary ID for match.
        >> id:str [None] = desired ID.
        >> gnspacc:boolean [True] = whether to try generating GnSpAcc formats for testing.
        >> uniprot:boolean [True] = whether to try making "ID (Acc)" combos for testing.
        << Returns True if ID found or False if not.
        '''
        if not id: return False
        if id in [self.info['ID'],self.info['AccNum'],self.list['Secondary ID'],self.shortName()]: return True
        fullacc = [self.info['AccNum']] + self.list['Secondary ID']
        if gnspacc:
            gnsp = '%s_%s' % (self.info['Gene'],self.info['SpecCode'])
            for acc in fullacc:
                if id == '%s__%s' % (gnsp,acc): return True
        if uniprot:
            for acc in fullacc:
                if id in ['%s (%s)' % (self.info['ID'],acc),'%s_%s' % (acc,self.info['SpecCode'])]: return True
        return False
#########################################################################################################################
    def fudgeFactor(self,posdict={},case=False,gaps=True,wildcards=True):   ### Returns the necessary fudge factor to match posdict
        '''
        Returns the necessary fudge factor to match posdict. Will raise error if it cannot do it!
        >> posdict:dict = {pos (1 to L):sequence}
        >> case:bool [False] = whether to use self.dict['Case'] information.
        >> gaps:bool [True] = whether to leave gaps in the sequence
        >> wildcards:bool [True] = whether to account for masked/wildcard Xs in posdict sequence
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not posdict: return 0
            maxf = max(rje.sortKeys(posdict)[0],self.aaLen() - rje.sortKeys(posdict)[-1])
            f = 0   # Fudge factor
            sequence = self.getSequence(case,gaps)
            ### ~ [2] Check fudging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while f < maxf:
                ## ~ [2a] Try forwards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fudged = True
                for pos in posdict:
                    if wildcards: pos_re = string.replace(posdict[pos],'X','\S')
                    else: pos_re = posdict[pos]
                    if not re.search('^%s' % pos_re,sequence[pos-1+f:]):
                        fudged = False
                        break
                if fudged: return f
                ## ~ [2b] Try backwards ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                fudged = True
                for pos in posdict:
                    if wildcards: pos_re = string.replace(posdict[pos],'X','\S')
                    else: pos_re = posdict[pos]
                    if not re.search('^%s' % pos_re,sequence[pos-1-f:]):
                        fudged = False
                        break
                if fudged: return -f
                f += 1
        except: self.log.errorLog('Major problem with %s fudgeFactor()' % self.shortName())
        raise ValueError
#########################################################################################################################
    ### <4> ### Sequence calculation methods                                                                            #
#########################################################################################################################
    def disorder(self,returnlist=False,reset=False): ### Adds disorder object and predicts disorder. See rje_disorder.py for details.
        '''Adds disorder object and predicts disorder. See rje_disorder.py for details.'''
        ### ~ [1] ~ Check for existing disorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        done = self.obj['Disorder'] and self.obj['Disorder'].list['ResidueDisorder'] and not reset
        ### ~ [2] ~ Add Object and calculate disorder, if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not done:
            self.obj['Disorder'] = rje_disorder.Disorder(self.log,self.cmd_list)
            done = self.obj['Disorder'].disorder(sequence=self.info['Sequence'],name=self.info['Name'])
        ### ~ [3] ~ Return list of disordered residues, or simply whether the disorder prediction worked ~~~~~~~~~~~~ ###
        if returnlist:
            if done: return self.obj['Disorder'].list['ResidueDisorder'][0:]
            else: return []
        else: return done
#########################################################################################################################
    def isDisordered(self,pos=-1,reset=False): ### Returns whether a particular residue is disordered
        '''
        Returns whether a particular residue is disordered.
        >> pos:int [-1] = position of residue (0->(L-1)). If < 0 will return True/False for status only.
        >> reset:bool [False] = Reset disorder list before assessing.
        '''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if reset or 'IsDisordered' not in self.list or not self.list['IsDisordered']:
            if not self.disorder(reset=reset):
                self.errorLog('Disorder prediction failed for %s' % (self.shortName()),printerror=False)
                raise ValueError
            disorder = [True] * self.aaLen()
            self.list['IsDisordered'] = [False] * self.aaLen()
            for region in self.obj['Disorder'].list['RegionDisorder']:
                self.list['IsDisordered'] = self.list['IsDisordered'][:region[0]-1] + disorder[region[0]-1:region[1]] + self.list['IsDisordered'][region[1]:]
        ### ~ [2] ~ Return relevant information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if pos < 0: return len(self.list['IsDisordered']) > 0
        return self.list['IsDisordered'][pos]
#########################################################################################################################
    def gappedDisorder(self,gap='prev'):   ### Returns gapped disorder list
        '''Returns gapped disorder list.'''
        ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        disorder = self.disorder(returnlist=True)
        gapdis = []
        ### ~ [1] ~ Add gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for r in self.getSequence():
            if r == '-':
                if gapdis and gap == 'prev': gapdis.append(gapdis[-1])
                elif gap == 'prev': gapdis.append(disorder[0])
                else: gapdis.append(gap)
            else: gapdis.append(disorder.pop(0))
        return gapdis
#########################################################################################################################
    def globProportion(self,absolute=False):     ### Returns proportion defined as globular
        '''Returns proportion defined as globular.'''
        self.disorder()
        return self.obj['Disorder'].globProportion(absolute)
#########################################################################################################################
    def aaFreq(self,aafreq={},newkeys=True):    ### Adds to aafreq dictionary (if given) and returns
        '''
        Adds to aafreq dictionary (if given) and returns.
        >> aafreq:dictionary of {aa:freq(count)}
        >> newkeys:boolean [True] = whether to add new AA keys if missing from aafreq
        << aafreq:new dictionary of values
        '''
        try: return aaFreq(self.info['Sequence'],aafreq,newkeys)
        except:
            self.log.errorLog('Problem with rje_sequence.aaFreq(%s)' % self.shortName(),quitchoice=True)
            return aafreq
#########################################################################################################################
    def sixFrameTranslation(self):  ### Translates DNA in all six reading frames into 'Translation' dictionary
        '''Translates DNA in all six reading frames into 'Translation' dictionary.'''
        try:self.dict['Translation'] = sixFrameTranslation(self.info['Sequence'][0:])
        except: self.log.errorLog('Problem with Sequence.sixFrameTranslation()')
#########################################################################################################################
    ### <5> ### Sequence masking methods                                                                                #
#########################################################################################################################
    def maskAA(self,maskaa,mask='X',log=True):    ### Masks given residues by type
        '''
        Adds disorder object with prediction and masks disorder.
        >> maskaa:list of AAs to be masked
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        << returns number of AAs masked
        '''
        try:### ~ [1] ~ Straight masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0
            for aa in maskaa:
                mx += self.info['Sequence'].count(aa)
                self.info['Sequence'] = self.info['Sequence'].replace(aa,mask)
            if log: self.printLog('#MASK','AA Mask (%s): %s %s added to %s' % (string.join(maskaa,';'),rje.integerString(mx),mask,self.shortName()))
            return mx
        except: self.errorLog('Problem masking AAs in %s' % (self.shortName())); return 0
#########################################################################################################################
    def maskDisorder(self,inverse=False,mask='X',log=True):   ### Adds disorder object with prediction and masks disorder
        '''
        Adds disorder object with prediction and masks disorder.
        >> inverse:bool [False] = Masks out predicted ordered regions (i.e. keeps disorder only)
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        '''
        try:### ~ [1] ~ Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deGap()    #!# Update to ignore gaps at some point? #!#
            if not self.disorder(): return self.errorLog('Disorder prediction failed for %s' % (self.shortName()),printerror=False)
            oldseq = self.info['Sequence'][0:]
            newseq = mask * len(oldseq)
            prex = oldseq.count(mask)
            if inverse: (newseq,oldseq) = (oldseq,newseq)
            ### ~ [2] ~ Mask relevant regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for region in self.obj['Disorder'].list['RegionDisorder']: oldseq = oldseq[:region[0]-1] + newseq[region[0]-1:region[1]] + oldseq[region[1]:]
            #x#self.deBug('%s::%s\n%s' % (self.obj['Disorder'].list['RegionDisorder'],self.obj['Disorder'].list['RegionFold'],oldseq))               
            ### ~ [3] ~ Update sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if maskx > 0:
                if inverse: mtxt = '%s masked %d ' % (self.shortName(),len(self.obj['Disorder'].list['RegionFold']))
                else: mtxt = '%s masked %d dis' % (self.shortName(),len(self.obj['Disorder'].list['RegionDisorder']))
                self.printLog('#MASK','%sordered regions. (%d %s added.)' % (mtxt,maskx,mask),screen=log)
        except: self.errorLog('Problem masking disorder from %s' % (self.shortName()))
#########################################################################################################################
    def maskLowComplexity(self,lowfreq=5,winsize=10,mask='X',log=True):     ### Masks low complexity regions of sequence
        '''
        Masks low complexity regions of sequence.
        >> lowfreq:int = Number of same aas in window size to mask
        >> winsize:int = Size of window to consider
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        '''
        try:### Setup ###
            self.deGap()
            if lowfreq > winsize or lowfreq < 3: return
            oldseq = self.info['Sequence'][0:]
            prex = oldseq.count(mask)

            ### Mask ###
            for r in range(len(oldseq)):
                a = self.info['Sequence'][r]
                if a == mask: continue
                x = 1
                for i in range(1,winsize):
                    if (r+i) >= len(oldseq): break
                    if self.info['Sequence'][r+i] == a: x += 1
                    if x >= lowfreq:
                        oldseq = oldseq[:r+1] + string.replace(oldseq[r+1:r+i],a,mask) + oldseq[r+i:]
                        #X#self.deBug('%s => %s' % (self.info['Sequence'][r+1:r+i],string.replace(oldseq[r+1:r+i],a,mask)))
                        break
                
            ### Update ###
            self.info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if maskx > 0: self.log.printLog('#MASK','Masked %s "low complexity" regions. (%d %s added.)' % (self.shortName(),maskx,mask),screen=log)
                    
        except: self.log.errorLog('Problem masking low complexity for %s' % self.shortName())
#########################################################################################################################
    def maskCase(self,case='lower',mask='X',log=True):   ### Masks sequence of given case
        '''
        Returns sequence.
        >> case:str ['lower'] = whether mask lower/upper case information.
        >> mask:str ['X'] = character to replace sequence with
        >> log:bool [True] = whether to log the amount of masking
        '''
        try:
            ### Setup ###
            if case.lower()[:1] == 'l': case = 'Lower'
            elif case.lower()[:1] == 'u': case = 'Upper'
            else: return False

            ### Mask ###
            prex = self.info['Sequence'].upper().count(mask.upper())
            sequence = self.info['Sequence'].upper()
            for (start,stop) in self.dict['Case'][case]: sequence = sequence[:start] + mask * (stop + 1 - start) + sequence[(stop+1):]
            self.info['Sequence'] = sequence

            ### Finish ###
            maskx = self.info['Sequence'].upper().count(mask.upper()) - prex
            if maskx > 0: self.log.printLog('#MASK','Masked %s %s case regions. (%d %s added.)' % (self.shortName(),case,maskx,mask),screen=log)
        except:
            self.log.errorLog('Problem masking case "%s" for %s' % (case,self.shortName()))
#########################################################################################################################
    def maskRegion(self,maskregion=[],inverse=False,mask='X',log=True): ### Masks region(s) of protein
        '''
        Masks region of protein.
        >> maskregion:list [] = Pairs of positions to (inclusively) mask
        >> inverse:bool [False] = Masks out predicted ordered regions (i.e. keeps disorder only)
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to print masking to log
        '''
        try:### ~ [1] ~ Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not maskregion: return 0
            maskregion = maskregion[0:]
            self.deGap()    #!# Update to ignore gaps at some point? #!#
            oldseq = self.info['Sequence'][0:]
            newseq = mask * len(oldseq)
            prex = oldseq.count(mask)
            if inverse: (newseq,oldseq) = (oldseq,newseq)
            ### ~ [2] ~ Mask relevant regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0
            while maskregion:
                start = maskregion.pop(0)
                try: end = maskregion.pop(0)
                except: end = -1
                if start < 0: start = self.aaLen() - start
                if end < 0: end = self.aaLen() - end
                end += 1
                oldseq = oldseq[:start] + newseq[start:end] + oldseq[end:]
                mx += 1
            ### ~ [3] ~ Update sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if maskx > 0: self.printLog('#MASK','%d region(s). (%d %s added.)' % (mx,maskx,mask),screen=log)
            return maskx
        except: self.errorLog('Problem with Region masking for %s' % (self.shortName()))
#########################################################################################################################
    def maskPosAA(self,maskdict={},mask='X',log=True):  ### Masks position-specific amino acids
        '''
        Masks position-specific amino acids. Returns and updates sequence.
        >> maskdict:dictionary of {pos:AAs} where pos is 1->L and AAs is a string of the AAs to mask
        >> mask:str ['X'] = character to replace sequence with
        >> log:bool [True] = whether to log the amount of masking
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = self.info['Sequence'].upper().count(mask.upper())
            sequence = self.info['Sequence'].upper()[0:]

            ### ~ [2] Mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pos in rje.sortKeys(maskdict):
                aas = maskdict[pos].upper()
                i = int(pos)
                if i > 0: i -= 1    # Can have backwards counting too and that's the same!
                try:
                    if aas in ['*','X','.'] or sequence[i] in rje.strList(aas): sequence = rje.strSub(sequence,i,i,mask)
                except: pass
            self.info['Sequence'] = sequence

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maskx = self.info['Sequence'].upper().count(mask.upper()) - prex
            if maskx > 0: self.log.printLog('#MASK','Masked %s position-specific AAs. (%d %s added.)' % (self.shortName(),maskx,mask),screen=log)
        except: self.log.errorLog('Problem masking positions for %s' % (self.shortName()))
        return self.info['Sequence']
#########################################################################################################################
## End of SECTION II: Sequence Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def specCodeFromName(name):  ### Returns the species code from a sequence name
    '''
    Returns the species code from a sequence name. (Stripped down version of Sequence.extractDetails())
    >> name:str = Sequence name
    << spcode:str = Species code
    '''
    seqname = string.replace(name,'||','|')  #!# Dodgy!
    accnum = seqname
    spcode = 'UNK'
    if rje.regExp(re_gnspacc,seqname):
        (gene,spcode,accnum) = rje.regExp(re_gnspacc,seqname)
    elif rje.regExp(re_gn_sp__acc,seqname):
        (gene,spcode,accnum) = rje.regExp(re_gn_sp__acc,seqname)
    elif rje.regExp(re_uniprot,seqname):
        (db,accnum,gene,spcode) = rje.regExp(re_uniprot,seqname)
    elif rje.regExp(re_uniprot2,seqname):
        (gene,spcode,accnum) = rje.regExp(re_uniprot2,seqname)
    elif rje.regExp(re_ipi,seqname):
        (accnum,_taxa) = rje.regExp(re_ipi,seqname)
        if _taxa in ipi_taxa.keys():
            spcode = ipi_taxa[_taxa]
        else:
            spcode = 'UNK'
    elif rje.regExp(re_ncbi,seqname):
        (id,db,accnum) = rje.regExp(re_ncbi,seqname)
        if rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',seqname):
            [accnum,gene,spcode] = rje.matchExp('\|sp\|(\S+)\|(\S+)_(\S+)',seqname)
        else:
            if rje.matchExp('\[(.+)\]\s*$',seqname):
                species = rje.matchExp('\[(.+)\]\s*$',seqname)[0]
                while re.search('\[(.+)$',species):
                    species = rje.matchExp('\[(.+)$',species)[0]
                spcode = getSpecCode(species)
    elif rje.regExp(re_ncbi2,seqname):
        (id,db,accnum) = rje.regExp(re_ncbi2,seqname)
        if rje.matchExp('\[(.+)\]\s*$',seqname):
            species = rje.matchExp('\[(.+)\]\s*$',seqname)[0]
            while re.search('\[(.+)$',species): species = rje.matchExp('\[(.+)$',species)[0]
            spcode = getSpecCode(species)
    elif re_genbank.match(seqname): (db,accnum) = rje.regExp(re_genbank,seqname)
    elif re_unigene.match(seqname): accnum = rje.regExp(re_unigene,seqname)[0]
    elif re_flybase.match(seqname):
        (acc,gene,spec) = rje.regExp(re_flybase,seqname)
        spcode = '%sRO%s' % (spec[0].upper(),spec[1:3].upper())
    elif re_ensemblpep.match(seqname):
        accnum = rje.regExp(re_ensemblpep,seqname)[0]
        spcode = 'UNK'
        if accnum.find('ENSAPMP') == 0: spcode = 'APIME'
        elif accnum.find('ENSBTAP') == 0: spcode = 'BOVIN'
        elif accnum.find('ENSGALP') == 0: spcode = 'CHICK'
        elif accnum.find('ENSPTRP') == 0: spcode = 'PANTR'
        elif accnum.find('ENSCINP') == 0: spcode = 'CIOIN'
        elif accnum.find('ENSCAV') == 0: spcode = 'CIOSA'
        elif accnum.find('ENSCAFP') == 0: spcode = 'CANFA'
        elif accnum.find('CG') == 0: spcode = 'DROME'
        elif accnum.find('ENSDNOP') == 0: spcode = 'DASNO'
        elif accnum.find('ENSETEP') == 0: spcode = 'ECHTE'
        elif accnum.find('SINFRUP') == 0: spcode = 'FUGRU'
        elif accnum.find('ENSGACP') == 0: spcode = 'GASAC'
        elif accnum.find('ENSDARP') == 0: spcode = 'BRARE'
        elif accnum.find('ENSP0') == 0: spcode = 'HUMAN'
        elif accnum.find('ENSLAFP') == 0: spcode = 'LOXAF' # Elephant
        elif accnum.find('ENSANGP') == 0: spcode = 'ANOGA'
        elif accnum.find('ENSMUSP') == 0: spcode = 'MOUSE'
        elif accnum.find('ENSOCUP') == 0: spcode = 'RABIT'  # Rabbit
        elif accnum.find('ENSRNOP') == 0: spcode = 'RAT'
        elif accnum.find('GSTENP') == 0: spcode = 'TETNG'
        elif accnum.find('ENSXETP') == 0: spcode = 'XENLA'
        elif accnum.find('AAEL') == 0: spcode = 'AEDAE'
        elif seqname.find('chromosome:SGD') > 0: spcode = 'YEAST'
        elif self.info['Name'].find('chromosome:CEL') > 0: spcode = 'CAEEL'
        elif self.info['Name'].find('ENSMMUP') == 0: spcode = 'MACMU'
        elif self.info['Name'].find('scaffold:MMUL') > 0: spcode = 'MACMU'
        elif self.info['Name'].find('ENSMODP') == 0: spcode = 'MONDO'
        elif self.info['Name'].find('scaffold:BROAD') > 0: spcode = 'MONDO'
    elif re_pipe.match(seqname): spcode = rje.regExp(re_pipe,seqname)[1]

    ## Special Databases ##
    if seqname.find('TP_') >= 0 and seqname.find('{Treponema pallidum Nichols}') >= 0: spcode = 'TREPA'
    if spcode == 'UNK' and (seqname.find('Emiliania huxleyi cDNA') >= 0 or accnum[:4] == 'EHUX'): spcode = 'EMIHU'
    if seqname.find('jgi|Thaps3|') == 0: spcode = 'THAPS'

    return spcode
#########################################################################################################################
def getSpecCode(species):   ### Returns spec_code for given species
    '''Returns spec_code for given species. This should be moved later and expanded to allow to read from speclist.'''
    spec_dic = {'Mus musculus':'MOUSE',
                'Homo sapiens':'HUMAN',
                'Gallus gallus':'CHICK',
                'Rattus norvegicus':'RAT',
                'Bos taurus':'BOVIN',
                'Danio rerio':'BRARE',
                'Saccharomyces cerevisiae':'YEAST'}
    if species in spec_dic.keys(): return spec_dic[species]
    if rje.matchExp('^(\S\S\S)\S*\s+(\S\S)',species):
        code = rje.matchExp('^(\S\S\S)\S*\s+(\S\S)',species.upper())
        return '%s%s' % (code[0],code[1])
    else: return string.replace(species,' ','').upper()[:10]
#########################################################################################################################
def eisenbergHydropathy(sequence,returnlist=False):  ### Returns the Eisenberg Hydropathy for the sequence
    '''
    Returns the Eisenberg Hydropathy for the sequence.
    >> sequence:str = AA sequence
    >> returnlist:boolean [False] = Returns list of hydropathies rather than total
    << hyd:float = Eisenberg Hydropathy for the sequence
    '''
    ### Setup ###
    eishyd = {'A':0.62, 'R':-2.53, 'N':-0.78, 'D':-0.9, 'C':0.29, 'Q':-0.85, 'E':-0.74, 'G':0.48, 'H':-0.4, 'I':1.38,
              'L':1.06, 'K':-1.5, 'M':0.64, 'F':1.19, 'P':0.12, 'S':-0.18, 'T':-0.05, 'W':0.81, 'Y':0.26, 'V':1.08,
              '-':0, 'X':0}
    eislist = []
    ### Calculate ###
    for aa in sequence: eislist.append(eishyd[aa])
    if returnlist: return eislist
    return sum(eislist)
#########################################################################################################################
def aaFreq(sequence,aafreq={},newkeys=True):    ### Adds to aafreq dictionary (if given) and returns
    '''
    Adds to aafreq dictionary (if given) and returns.
    >> sequence:str = Sequence for AAFreq calculation
    >> aafreq:dictionary of {aa:freq(count)}
    >> newkeys:boolean [True] = whether to add new AA keys if missing from aafreq
    << aafreq:new dictionary of values
    '''
    for aa in sequence:
        if aafreq.has_key(aa): aafreq[aa] += 1
        elif newkeys: aafreq[aa] = 1.0
    return aafreq
#########################################################################################################################
def trypDigest(sequence=''):    ### Returns trypsin digestion of sequence
    '''Returns trypsin digestion of sequence.'''
    digest = string.join(string.split(sequence.upper(),'K'),'K ')
    digest = string.join(string.split(digest,'R'),'R ')
    digest = string.split(digest)
    return digest
#########################################################################################################################
def MWt(sequence=''):   ### Returns Molecular Weight of Sequence
    '''Returns Molecular Weight of Sequence.'''
    mwtdic = {'A':89, 'V':117, 'L':131, 'I':131, 'P':115,
              'F':165, 'W':204, 'M':149, 'G':75, 'S':105,
              'T':119, 'C':121, 'Y':181, 'N':132, 'Q':146,
              'D':133, 'E':147, 'K':146, 'R':174, 'H':155,
              'X':136.75}
    _mwt = 0.0
    for aa in sequence.upper():
        if aa in mwtdic.keys():
            _mwt += mwtdic[aa]
            _mwt += 18  # H2O
    if _mwt > 18:
        _mwt -= 18
    return _mwt
#########################################################################################################################
def chargeDict(sequence,callobj=None):   ### Performs absolute, net and charge balance calculations 
    '''
    Performs absolute, net and charge balance calculations.
    >> sequence:str = sequence to calculate stats on.
    >> callobj:Object = calling object for error messages etc.
    << dictionary of {Type:Value)
    '''
    try:
        charge = []
        for a in sequence: # Motif
            if a in ['K','R']:
                charge.append(1)
            elif a in ['D','E']:
                charge.append(-1)
            else:
                charge.append(0)
        return {'AbsChg':(charge.count(1) + charge.count(-1)),
                'NetChg':sum(charge),
                'BalChg':(sum(charge[:int(len(charge)/2)]) - sum(charge[-int(len(charge)/2):]))}
    except:
        if callobj:
            callobj.log.errorLog('Error in rje_sequence.chargeDict()')
            return {}
        raise
#########################################################################################################################
def peptideDetails(sequence,callobj=None):  ### Returns OK if sequence alright, or warning if bad aa combos
    '''
    Returns OK if sequence alright, or warning if bad aa combos.
    >> sequence:str = sequence to calculate stats on.
    >> callobj:Object = calling object for error messages etc.
    << string of peptide assessment
    '''
    try:
        bad_combo = []
        for combo in ['DP','DC','DG','NG','NS','PPP']:
            if string.strip(sequence.upper(),'-').find(combo) >= 0: bad_combo.append(combo)
        if bad_combo: return 'Warning: %s' % string.join(bad_combo,', ')
        else: return 'OK'
    except:
        if callobj:
            callobj.log.errorLog('Error in rje_sequence.peptideDetails()')
            return {}
        raise
#########################################################################################################################
def mapGaps(inseq,gapseq,callobj=None):  ### Returns inseq with gaps inserted as found in gapseq
    '''Returns inseq with gaps inserted as found in gapseq.'''
    try:
        newseq = ''
        while len(newseq) < len(gapseq):
            if gapseq[len(newseq)] == '-': newseq += '-'
            elif inseq:
                newseq += inseq[:1]
                inseq = inseq[1:]
            else:
                if callobj: callobj.log.errorLog('Sequence length mismatch during mapGaps()',printerror=False)
                while len(newseq) < len(gapseq): newseq += '-'
        return newseq
    except:
        print newseq, inseq
        print gapseq
        return newseq
#########################################################################################################################
def dna2prot(dnaseq):   ### Returns a protein sequence for a given DNA sequence
    '''Returns a protein sequence for a given DNA sequence.'''
    prot = ''
    while len(dnaseq)>2:
        codon = dnaseq[:3]
        dnaseq = dnaseq[3:]
        try: prot += genetic_code[string.replace(codon.upper(),'T','U')]
        except: prot += 'X'
    return prot
#########################################################################################################################
def reverseComplement(dnaseq,rna=False):  ### Returns the reverse complement of the DNA sequence given (upper case)
    '''Returns the reverse complement of the DNA sequence given.'''
    revcomp = ''
    repdict = {'G':'C','C':'G','T':'A','A':'T'}
    if rna: repdict['A'] = 'U'
    for aa in string.replace(dnaseq.upper(),'U','T'):   # Convert RNA to DNA
        try: revcomp = '%s%s' % (repdict[aa],revcomp)    
        except: revcomp = '%s%s' % ('N',revcomp)    
    return revcomp 
#########################################################################################################################
def sixFrameTranslation(dnaseq):  ### Translates DNA in all six reading frames into dictionary
    '''Translates DNA in all six reading frames into dictionary.'''
    sixrf = {}
    for rf in range(3): sixrf[rf+1] = dna2prot(dnaseq[rf:])
    dnaseq = reverseComplement(dnaseq)
    for rf in range(3): sixrf[-(rf+1)] = dna2prot(dnaseq[rf:])
    return sixrf
#########################################################################################################################
def threeFrameTranslation(dnaseq,minpoly=0):  ### Translates DNA in all three reading frames into dictionary
    '''Translates DNA in all six reading frames into dictionary.'''
    a = 1
    while minpoly >0 and dnaseq[-a].upper() == 'A': a += 1
    a -= 1
    if a >= minpoly and minpoly > 0: dnaseq = dnaseq[:-a]
    rf3 = {}
    for rf in range(3): rf3[rf+1] = dna2prot(dnaseq[rf:])
    return rf3
#########################################################################################################################
def bestORF(protseq,startm=False,nonx=True):  ### Returns longest ORF (or first of equal length)
    '''
    Returns longest ORF (or first of equal length).
    >> sequence:str = Protein sequence, i.e. translation of DNA/RNA.
    >> startm:bool [False] = Whether the ORF should start with a methionine
    >> nonx:bool [True] = Whether ORFs should only be assessed in terms of their non-X content.
    '''
    (bestorf,bestlen) = ('',0)
    for orf in string.split(protseq,'*'):
        if startm and orf.count('M') < 1: continue
        if startm and orf.find('M') >= 0: orf = orf[orf.find('M'):]
        orflen = len(orf)
        if nonx: orflen -= string.count(orf.upper(),'X')
        if orflen > bestlen: (bestorf,bestlen) = (orf,orflen)
    return bestorf
#########################################################################################################################
def estTranslation(dnaseq,minpoly=10):  ### Returns translations of EST into protein reading frames as dictionary
    '''
    Returns translation of EST into protein reading frames.
    >> dnaseq:str = DNA sequence to translate
    >> minpoly:str = Min. length of poly-A or poly-T to be recognised and removed.
    '''
    dnaseq = dnaseq.upper()
    t = 0       # Look for run of 5' Ts (i.e. rev complement of poly-A tail
    while dnaseq[t] == 'T': t += 1
    a = 1
    while dnaseq[-a] == 'A': a += 1
    a -= 1
    poly = ''
    if a >= minpoly and minpoly >= 0:
        dnaseq = dnaseq[:-a]
        poly = 'A'
    elif t >= minpoly and minpoly >= 0:
        dnaseq = dnaseq[t:]
        poly = 'T'
    sixrf = sixFrameTranslation(dnaseq)
    if poly in ['A','T']:
        for rf in rje.sortKeys(sixrf):
            if rf > 0  and poly == 'T': sixrf.pop(rf)
            elif rf < 0  and poly == 'A': sixrf.pop(rf)
            else: sixrf[rf] = '%s*' % string.join(string.split(sixrf[rf],'*')[:-1],'*')
    return sixrf
#########################################################################################################################
def estTrunc(dnaseq,minpoly=10):  ### Returns truncation of EST along with poly-AT as tuple (5',seq,3')
    '''
    Returns translation of EST into protein reading frames.
    >> dnaseq:str = DNA sequence to translate
    >> minpoly:str = Min. length of poly-A or poly-T to be recognised and removed.
    '''
    dnaseq = dnaseq.upper()
    t = 0       # Look for run of 5' Ts (i.e. rev complement of poly-A tail
    while dnaseq[t] == 'T': t += 1
    a = 1
    while dnaseq[-a] == 'A': a += 1
    a -= 1
    poly = ''
    if a >= minpoly and minpoly >= 0: return ('',dnaseq[:-a],dnaseq[-a:])
    elif t >= minpoly and minpoly >= 0: return (dnaseq[:t],dnaseq[t:],'')
    return ('',dnaseq,'')
#########################################################################################################################
def caseDict(sequence): ### Return dictionary of {'Upper':[(start,end)],'Lower':[(start,end)]}
    '''Return dictionary of {'Upper':[(start,end)],'Lower':[(start,end)]}.'''
    case = {'Upper':[],'Lower':[]}      # Dictionary
    r = 0                               # Residue position
    for block in re.findall(re.compile('[A-Z\-\*]+|[a-z\-\*]+'),sequence):
        ckey = 'Lower'
        if block == block.upper(): ckey = 'Upper'
        case[ckey].append((r,r+len(block)-1))
        r += len(block)
    return case
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: print 'Not for standalone running!'
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV
#########################################################################################################################
