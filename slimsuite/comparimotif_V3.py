#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       comparimotif
Description:  Motif vs Motif Comparison Software
Version:      3.5
Last Edit:    04/08/09
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    CompariMotif is a piece of software with a single objective: to take two lists of protein motifs and compare them to
    each other, identifying which motifs have some degree of overlap, and identifying the relationships between those
    motifs. It can be used to compare a list of motifs with themselves, their reversed selves, or a list of previously
    published motifs, for example (e.g. ELM (http://elm.eu.org/)). CompariMotif outputs a table of all pairs of matching
    motifs, along with their degree of similarity (information content) and their relationship to each other.

    The best match is used to define the relationship between the two motifs. These relationships are comprised of the
    following keywords:

    Match type keywords identify the type of relationship seen:
    * Exact = all the matches in the two motifs are precise
    * Variant = focal motif contains only exact matches and subvariants of degenerate positions in the other motif
    * Degenerate = the focal motif contains only exact matches and degenerate versions of positions in the other motif
    * Complex = some positions in the focal motif are degenerate versions of positions in the compared motif, while 
    others are subvariants of degenerate positions.

    Match length keywords identify the length relationships of the two motifs:
    * Match = both motifs are the same length and match across their entire length
    * Parent = the focal motif is longer and entirely contains the compared motif
    * Subsequence = the focal motif is shorter and entirely contained within the compared motif
    * Overlap = neither motif is entirely contained within the other

    This gives sixteen possible classifications for each motif's relationship to the compared motif.

Input:
    CompariMotif can take input in a number of formats. The preferred format is SLiMSearch format, while is a single line
    motif format: 'Name Sequence #Comments' (Comments are optional and ignored). Alternative inputs include SLiMDisc and
    Slim Pickings output, raw lists of motifs, and fasta format.

Output:
    The main output for CompariMotif is delimited text file containing the following fields:
    * File1     = Name of motifs file (if outstyle=multi)
    * File2     = Name of searchdb file (if outstyle=multi)
    * Name1     = Name of motif from motif file 1 
    * Name2     = Name of motif from motif file 2 
    * Motif1    = Motif (pattern) from motif file 1
    * Motif2    = Motif (pattern) from motif file 2
    * Sim1      = Description of motif1's relationship to motif2
    * Sim2      = Description of motif2's relationship to motif1
    * Match     = Text summary of matched region
    * MatchPos  = Number of matched positions between motif1 and motif2 (>= mishare=X)
    * MatchIC   = Information content of matched positions
    * NormIC    = MatchA as a proportion of the maximum possible MatchA
    * Score     = Heuristic score (MatchPos x NormIC) for ranking motif matches
    * Info1     = Ambiguity score of motif1 
    * Info2     = Ambiguity score of motif2 
    * Desc1     = Description of motif1 (if motdesc = 1 or 3)
    * Desc2     = Description of motif2 (if motdesc = 2 or 3)

    With the exception of the file names, which are only output if outstyle=multi, the above is the output for the
    default "normal" output style. If outstyle=single then only statistics for motif2 (the searchdb motif) are output
    as this is designed for searches using a single motif against a motif database. If outstyle=normalsplit or
    outstyle=multisplit then motif1 information is grouped together, followed by motif2 information, followed by the
    match statistics. More information can be found in the CompariMotif manual.

Commandline:
    ## Basic Input Parameters ##
    motifs=FILE     : File of input motifs/peptides [None]
    searchdb=FILE   : (Optional) second motif file to compare. Will compare to self if none given. [None]
    dna=T/F         : Whether motifs should be considered as DNA motifs [False]

    ## Basic Output Parameters ##
    resfile=FILE    : Name of results file, FILE.compare.txt. [motifsFILE-searchdbFILE.compare.txt]
    motinfo=FILE    : Filename for output of motif summary table (if desired) [None]
    motific=T/F     : Output Information Content for motifs [False]

    ## Motif Comparison Parameters ##
    minshare=X      : Min. number of non-wildcard positions for motifs to share [2]
    normcut=X       : Min. normalised MatchIC for motif match [0.5]
    matchfix=X      : If >0 must exactly match *all* fixed positions in the motifs from:  [0]
                        - 1: input (motifs=FILE) motifs
                        - 2: searchdb motifs
                        - 3: *both* input and searchdb motifs
    ambcut=X        : Max number of choices in ambiguous position before replaced with wildcard (0=use all) [10]

    ## Advanced Motif Input Parameters ##
    minic=X         : Min information content for a motif (1 fixed position = 1.0) [2.0]
    minfix=X        : Min number of fixed positions for a motif to contain [0]
    minpep=X        : Min number of defined positions in a motif [2]
    trimx=T/F       : Trims Xs from the ends of a motif [False]
    nrmotif=T/F     : Whether to remove redundancy in input motifs [False]
    reverse=T/F     : Reverse the input motifs. [False]
                        - If no searchdb given, these will be searched against the "forward" motifs.
    mismatches=X    : <= X mismatches of positions can be tolerated [0]
    aafreq=FILE     : Use FILE to replace uniform AAFreqs (FILE can be sequences or aafreq) [None]    

    ## Advanced Motif Output Parameters ##
    xgmml=T/F       : Whether to output XGMML format results [True]
    pickle=T/F      : Whether to load/save pickle following motif loading/filtering [False]
    motdesc=X       : Sets which motifs have description outputs (0-3 as matchfix option) [3]
    outstyle=X      : Sets the output style for the resfile [normal]
                        - normal = all standard stats are output
                        - multi = designed for multiple appended runs. File names are also output
                        - single = designed for searches of a single motif vs a database. Only motif2 stats are output
                        - reduced = motifs do not have names or descriptions
                        - normalsplit/multisplit = as normal/multi but stats are grouped by motif rather than by type

Uses general modules: os, string, sys, time
Uses RJE modules: rje, rje_menu, rje_slim, rje_slimlist, rje_xgmml, rje_zen
Other modules needed: rje_aaprop, rje_dismatrix, rje_seq, rje_blast, rje_pam, rje_sequence, rje_uniprot
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, sets, string, sys, time
#########################################################################################################################
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_menu, rje_seq, rje_slim, rje_slimlist, rje_xgmml, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Full working version with menu
    # 1.1 - Added extra output options
    # 2.0 - Reworked for functionality with MotifList instead of PRESTO and using own methods.
    # 2.1 - Minor bug fixing and tidying. Removed matchic=F option. Added score and normcut=X.
    # 3.0 - Replaced rje_motif* modules with rje_slim* modules and improved handling of termini.
    # 3.1 - Added XGMML output.
    # 3.2 - Added mismatches=X option. (NB. mismatch=X is used in SLiMSearch.)
    # 3.3 - Added "Match" column, summarising matches
    # 3.4 - Added a DNA option and AA frequencies.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    [Y] : Update the docstring to account for new module use.
    [ ] : Add weighted IC and use of AAfreq files
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('CompariMotif', '3.5', 'August 2009', '2007')
    description = 'Motif vs Motif Comparison'
    author = 'Dr Richard J. Edwards.'
    comments = ['Cite: Edwards RJ, Davey NE & Shields DC (2008). Bioinformatics 24(10):1307-9.']
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show rje_slimlist options?'): out.verbose(-1,4,text=rje_slimlist.__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: CompariMotif Class                                                                                      #                                                                                                          #
#########################################################################################################################
class CompariMotif(rje.RJE_Object):     
    '''
    CompariMotif Motif-Motif Comparison Class. Author: Rich Edwards (2005).

    Info:str
    - AAFreq = Use FILE to replace uniform AAFreqs (FILE can be sequences or aafreq) [None]    
    - MotInfo = Filename for output of motif summary table (if desired) [None]
    - OutStyle = Output Style for comparison results [normal]
    - SearchDB = Second file of motifs to compare []
        
    Opt:boolean
    - DNA = Whether motifs should be considered as DNA motifs [False]
    - MotifIC  Output Information Content for motifs [False]
    - Pickle = Whether to load/save pickle following motif loading/filtering [False]
    - XGMML = Whether to output XGMML format results [True]

    Stat:numeric
    - NormCut = Min. normalised MatchIC for motif match [0.5]
    - MatchFix = If >0 must match *all* fixed positions in the motifs from:  [0]
        - 1: input (motifs=FILE) motifs
        - 2: searchdb motifs
        - 3: *both* input and searchdb motifs
    - MinShare = Min. number of non-wildcard positions for motifs to share [2]
    - Mismatches = <= X mismatches of positions can be tolerated [0]
    - MotDesc = Sets which motifs have description outputs (0-3 as matchfix option) [3]

    List:list
    - Matches = List of match dictionaries used for Cytoscape output

    Dict:dictionary
    - AAFreq = Dictionary of AA Frequencies (uniform by default)

    Obj:RJE_Objects
    - SlimList = SLiMList object, which controls most command-line parameters!
    - SearchDB = SLiMList object containing compared dictionary
    '''
    def alphabet(self): return self.obj['SlimList'].list['Alphabet']
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AAFreq','MotInfo','OutStyle','SearchDB']
        self.optlist = ['MotifIC','DNA','XGMML','Pickle']
        self.statlist = ['MatchFix','NormCut','MinShare','Mismatches']
        self.listlist = ['Matches']
        self.dictlist = ['AAFreq']
        self.objlist = ['SlimList','SearchDB']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### SlimList ###
        mdefaults = ['nrmotif=F','ambcut=10']                   # Reset SlimList defaults
        self.cmd_list = mdefaults + self.cmd_list + ['compare=T']    # Compare used by Motif object!
        self.setStat({'MatchFix':0,'MotDesc':3,'NormCut':0.5,'MinShare':2,'Mismatches':0})
        self.setOpt({'XGMML':True})
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,['pickle=F']+self.cmd_list)
        if self.obj['SlimList'].info['Motifs'].lower() in ['','none'] and self.obj['SlimList'].stat['Interactive'] == 0:
            self.obj['SlimList'].stat['Interactive'] = 1
        self.setInfo(self.obj['SlimList'].info)
        self.setOpt(self.obj['SlimList'].opt)
        self.setStat(self.obj['SlimList'].stat)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Specific Commands ###
                self._cmdReadList(cmd,'info',['OutStyle'])
                self._cmdReadList(cmd,'opt',['MotifIC','DNA','XGMML'])
                self._cmdReadList(cmd,'file',['MotInfo','SearchDB','AAFreq'])
                self._cmdRead(cmd,type='info',att='Basefile',arg='resfile')
                self._cmdReadList(cmd,'int',['MatchFix','MotDesc','MinShare','Mismatches'])
                self._cmdReadList(cmd,'stat',['NormCut'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        ### Special update of self.obj['SlimList'] for Menu ###
        self.obj['SlimList'].setInfo(self.info)
        self.obj['SlimList'].setStat(self.stat)
        self.obj['SlimList'].setOpt(self.opt)
#########################################################################################################################
    ### <2> ### Main Run Method                                                                                         #
#########################################################################################################################
    def run(self):         ### Main run method. Sets up parameters and calls comparison method.
        '''Main run method. Sets up parameters and calls comparison method. Returns XGMML.'''
        try:
            ### Setup ###
            self.setupRun()

            ### Menu ###
            if self.obj['SlimList'].stat['Interactive'] >= 1:
                menuhead = 'CompariMotif Main Menu'
                # (code,desc,type,key)      #!# May have to split this into SlimList and Self menus #!#
                menulist = [('','### INPUT OPTIONS ###','',''),
                            ('M','Motif file','info','Motifs'),
                            ('D','SearchDB file','info','SearchDB'),
                            ('R','Reverse input motifs','opt','Reverse'),
                            ('N','DNA input motifs','opt','DNA'),
                            ('','\n### OUTPUT OPTIONS ###','',''),
                            ('O','Name of results file','info','Resfile'),
                            ('T','Summary motif table','info','MotInfo'),
                            ('I','Output motif Information Content','opt','MotifIC'),
                            ('','\n### COMPARISON OPTIONS ###','',''),
                            ('S','Min. no. positions to share','stat','MinShare'),
                            ('C','Normalised IC cut-off','stat','NormCut'),
                            ('F','Set MatchFix parameter (0/1/2/3)','stat','MatchFix'),
                            ('A','Max no. choices in ambiguous position','stat','AmbCut'),
                            ('','\n### ADVANCED INPUT OPTIONS ###','',''),
                            ('MP','Min. non-wildcard length of motif/peptide','stat','MinPos'),
                            ('MF','Min. no. fixed positions','stat','MinFix'),
                            ('TX','Trims wildcards from the ends of motif','opt','TrimX'),
                            ('MM','No. mismatches to tolerate','int','Mismatches'),
                            ('AA','AA Frequency file','info','AAFreq'),
                            ('','\n### ADVANCED OUTPUT OPTIONS ###','',''),
                            ('OS','Output style','info','OutStyle'),
                            ('MD','Motif Description output','int','MotDesc'),
                            ('AP','Append output','opt','Append'),
                            ('','','return',True),
                            ('X','Exit menu and proceed.','return',True)]
                rje_menu.menu(self.obj['SlimList'],menuhead,menulist,choicetext='Please select (Blank to proceed):',changecase=True)
                for opt in ['MotifIC']:
                    self.obj['SlimList'].cmd_list.append('%s=%s' % (opt.lower(),self.obj['SlimList'].opt[opt]))        
                self.obj['SlimList'].cmd_list.append('ambcut=%d' % self.obj['SlimList'].stat['AmbCut'])
            self.setupRun(postmenu=True)

            ### Perform Motif Comparisons ###
            self.motifCompare()
            return self.compareXGMML()
        except KeyboardInterrupt: raise
        except:
            self.log.errorLog('Error in CompariMotif.run()')
            raise
#########################################################################################################################
    def setupRun(self,postmenu=False):     ### Sets up runtime parameters
        '''Sets up runtime parameters.'''
        try:
            ### ~ [1] ~ Setup Input Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motifs = self.obj['SlimList']
            if postmenu:
                pickledmotifs = motifs.unpickleMe()
                if pickledmotifs: motifs = pickledmotifs
                else: motifs.loadMotifs()
                if not motifs.slimNum(): sys.exit()
                if pickledmotifs and self.opt['Pickle']: motifs.pickleMe() 
            if motifs.info['SearchDB'].lower() in ['none','']: motifs.info['SearchDB'] = motifs.info['Motifs']
            
            ### ~ [2] ~ Setup AA/DNA Frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aafile = self.info['AAFreq']
            if postmenu and aafile.lower() not in ['','none']:
                seqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F','accnr=F','seqnr=F'])
                if os.path.exists(aafile):
                    try:
                        if open(aafile,'r').read()[:1] == '>': self.dict['AAFreq'] = seqlist.aaFreq(fromfile=aafile,alphabet=self.alphabet(),total=False)    # Returns proportions
                        else: self.dict['AAFreq'] = seqlist.aaFreq(loadfile=aafile,alphabet=self.alphabet(),total=False)    # Returns proportions
                    except: self.log.errorLog('Cannot use AAFreq File')
                else: self.log.errorLog('Cannot find AA Frequency / Sequence file "%s"' % self.info['AAFreq'],printerror=False)
            if postmenu and not self.dict['AAFreq']:
                self.dict['AAFreq'] = {}
                for aa in self.alphabet(): self.dict['AAFreq'][aa] = 1.0 / len(self.alphabet())

            ### ~ [3] ~ Setup Default Result Filenames ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if motifs.info['Basefile'].lower() in ['','none']:
                mname = rje.baseFile(motifs.info['Motifs'],strip_path=True,extlist=['txt','motif','motifs','fas'])
                sname = rje.baseFile(motifs.info['SearchDB'],strip_path=True,extlist=['txt','motif','motifs','fas'])
                if 'none' in [mname.lower(),sname.lower()]: motifs.info['Basefile'] = 'comparimotif'
                else: motifs.info['Basefile'] = '%s-%s' % (mname,sname)
                if string.count(motifs.info['Basefile'],',') > 0: # Lists of motifs!
                    if motifs.info['SearchDB'] == motifs.info['Motifs']: motifs.info['Basefile'] = 'compari-self'
                    else: motifs.info['Basefile'] = 'comparimotif'
            if 'Resfile' in motifs.info and motifs.info['Resfile'].lower() not in ['','none']: return
            motifs.info['Resfile'] = '%s.compare.%s' % (motifs.info['Basefile'],rje.delimitExt(rje.getDelimit(self.cmd_list,'\t')))
        except:
            self.log.errorLog('Error in setupRun()')
            raise
#########################################################################################################################
    ### <3> ### Motif-Motif Comparisons                                                                                 #
#########################################################################################################################
    def motifCompare(self):     ### Compares two lists of motifs against each other
        '''Compares two lists of motifs against each other'''
        try:
            ### Setup Motif Lists and Output Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Parameter settings # 
            self.info = self.obj['SlimList'].info
            self.stat = self.obj['SlimList'].stat
            self.opt = self.obj['SlimList'].opt
            motifs = self.obj['SlimList']
            compdb = self.obj['SearchDB'] = rje_slimlist.SLiMList(self.log,self.obj['SlimList'].cmd_list+['reverse=F'])
            compdb.loadMotifs(file=self.info['SearchDB'])
            ## Headers ##
            #hlist = [['MatchPos','MatchIC','NormIC','Score','Score2'], ['MotifFile','Name1','Motif1','Sim1'], ['SearchDB','Name2','Motif2','Sim2']]
            hlist = [['Match','MatchPos','MatchIC','NormIC','Score'], ['MotifFile','Name1','Motif1','Sim1'], ['SearchDB','Name2','Motif2','Sim2']]
            for x in [1,2]:
                if self.opt['MotifIC']: hlist[x].append('Info%d' % x)
                if self.stat['MotDesc'] in [x,3]: hlist[x].append('Desc%d' % x)
            if self.info['OutStyle'].find('reduced') == 0: hlist = [hlist[0][0:],hlist[1][2:4],hlist[2][2:4]]
            if self.info['OutStyle'] == 'multisplit': headers = hlist[1] + hlist[2] + hlist[0]
            if self.info['OutStyle'] == 'normalsplit': headers = hlist[1][1:] + hlist[2][1:] + hlist[0]
            if self.info['OutStyle'] == 'reducedsplit': headers = hlist[1][0:] + hlist[2][0:] + hlist[0]
            elif self.info['OutStyle'] == 'single': headers = hlist[2] + hlist[0]
            else:
                if self.info['OutStyle'].lower() in ['normal']: i = 1
                else: i = 0
                headers = []
                while i < len(hlist[1]) or i < len(hlist[2]):
                    if i < len(hlist[1]): headers.append(hlist[1][i])
                    if i < len(hlist[2]): headers.append(hlist[2][i])
                    if headers[-1] == 'Sim2': headers += hlist[0]
                    i += 1
            ## Output File ##
            delimit = rje.getDelimit(self.cmd_list,default='\t')
            resfile = self.info['Resfile']
            rje.delimitedFileOutput(self,resfile,headers,delimit,rje_backup=True)
            minshare = self.stat['MinShare']
            self.list['Matches'] = []

            ### Setup vardict {motif:varlist of lists} and var_ic {element:ic} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            vardict = {}    ### Dictionary of motif:list of variants = lists of elements
            var_ic = {}     ### Dictionary of motif element and IC of said element
            compx = motifs.motifNum() + compdb.motifNum()
            cx = 0.0
            for motif in motifs.motifs() + compdb.motifs():
                cx += 100.0
                self.log.printLog('\r#COMP','Varlist setup: %.2f%%' % (cx/compx),log=False,newline=False)
                motif.weightedIC(self.dict['AAFreq'])
                motif.searchDict()
                vardict[motif] = []
                for searchvar in motif.dict['Search'][0][0:]:
                    addvar = []
                    while searchvar:
                        if searchvar[0] == '[':
                            el = rje.strSort(searchvar[1:searchvar.find(']')])
                            searchvar = searchvar[searchvar.find(']')+1:]
                        else:
                            el = searchvar[:1]
                            searchvar = searchvar[1:]
                        if not var_ic.has_key(el):
                            var_ic[el] = rje_slim.weightedElementIC(el,self.dict['AAFreq'])    ### Can add aafreq and wild_pen etc. here if desired.
                        addvar.append(el)
                    vardict[motif].append(addvar[0:])
                #X#self.deBug('%s: %s' % (motif.info['Motifs'],vardict[motif]))
            if motifs.info['MotInfo'].lower() not in ['','none']: motifs.motifInfo()    # Do after IC adjustment

            ### Compare Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            compx = len(motifs.motifs()) * len(compdb.motifs())
            (cx,mx) = (0.0,0)
            for motif1 in motifs.motifs():
                seq1 = motif1.pattern()     # string.replace(string.join(motif1.list['PRESTO'],''),'X','.')
                for motif2 in compdb.motifs():
                    cx += 100.0
                    self.log.printLog('\r#COMP','Comparing %s & %s: %.2f%%' % (hlist[1][0],hlist[2][0],cx/compx),log=False,newline=False)
                    seq2 = motif2.pattern()     # string.replace(string.join(motif2.list['PRESTO'],''),'X','.')
                    match = []      # Give value to skip later stages
                    outdict = {'MotifFile':rje.baseFile(self.info['Name'],strip_path=True),#?#'SLiM1':motif1,'SLiM2':motif2,
                               'SearchDB':rje.baseFile(self.info['SearchDB'],strip_path=True),
                               'Name1':motif1.info['Name'], 'Name2':motif2.info['Name'],'Motif1':motif1.info['Sequence'],'Motif2':motif2.info['Sequence'],
                               'Desc1':motif1.info['Description'], 'Desc2':motif2.info['Description'],
                               'Info1':'%.2f' % motif1.stat['IC'], 'Info2':'%.2f' % motif2.stat['IC']}

                    ## 1. Straight comparison of sequences ##
                    if seq1 == seq2:
                        match = ['Exact Match','Exact Match',motif1.slimPos(),motif1.stat['IC'],motif1.pattern()]
                        if self.info['Name'] == self.info['SearchDB'] and motif1.info['Name'] == motif2.info['Name']: continue  # Self!
                        elif motif1.info['Name'] == motif2.info['Name']: match = ['Duplicate','Duplicate',motif1.slimPos(),motif1.stat['IC']]
                    elif self.speedskip(seq1,seq2): continue    # No chance of match!

                    ## 2. Precise subsequences (motif.list['PRESTO']) ##
                    elif motif1.slim(wings=True).find(motif2.slim(wings=True)) >= 0:
                        match = ['Exact Parent','Exact Subsequence',motif2.slimPos(),motif2.stat['IC'],rje_slim.patternFromCode(motif2.slim())]
                    elif motif2.slim(wings=True).find(motif1.slim(wings=True)) >= 0:
                        match = ['Exact Subsequence','Exact Parent',motif1.slimPos(),motif1.stat['IC'],rje_slim.patternFromCode(motif1.slim())]
                
                    ## 3. Sliding comparison of vardict[motif] ##
                    if not match:
                        comtype = {'Good':['match', 'ambmatch', 'ambdegen', 'ambvar', 'degen', 'variant', 'overlap'],
                                   'Bad':['badamb', 'ambmismatch', 'mismatch'],
                                   'OK':['wild', 'wildvar','wilddegen'],
                                   'Fix':['match'],
                                   'Degen':['ambdegen','wilddegen','degen','overlap'],
                                   'Variant':['ambvar','variant','wildvar']} 
                        ## Setup list of possible results ##
                        vmatch = ()     # Tuple of best match: (ic,assessment1,assessment2)
                        for v1 in vardict[motif1]:      # Variant of motif1
                            if len(v1) < minshare: continue
                            for v2 in vardict[motif2]:  # Variant of motif2
                                if len(v2) < minshare: continue
                                #x#self.deBug('%s v %s ...' % (v1,v2))
                                ## Goes through every sliding comparison of at least minshare in length ##
                                for x in range(-(len(v1) - minshare),len(v2) - minshare + 1):    # This is how far v1 must stick out from v2
                                    if x < 0:
                                        c1 = v1[-x:-x+len(v2)]      # Overlapping elements of motif1
                                        c2 = v2[:len(v1)+x]         # Overlapping elements of motif2
                                    else:
                                        c1 = v1[:len(v2)-x]         # Overlapping elements of motif1
                                        c2 = v2[x:x+len(v1)]        # Overlapping elements of motif2
                                    m1 = []     # Types of match motif1 vs motif2
                                    m2 = []     # Types of match motif2 vs motif1
                                    ic = 0.0
                                    #x#self.deBug('... %s vs %s' % (c1,c2))
                                    for i in range(len(c1)):
                                        ## Determine relationship of each position ##
                                        pos_relationship = self._positionRelationship(c1[i],c2[i])
                                        #x#self.deBug('... %s vs %s = %s' % (c1[i],c2[i],pos_relationship))
                                        #!# Could assess matchfix here for runtime speed (but more coding mess) #!#
                                        m1.append(pos_relationship[0])
                                        m2.append(pos_relationship[1])
                                        # Good = match, ambmatch, ambdegen, ambvar, degen, variant, overlap
                                        # Bad  = badamb, ambmismatch, mismatch
                                        # OK   = wild, wildvar, 
                                        ## IC of match ##
                                        if m1[-1] not in comtype['Bad']:        # Add lower IC if match not bad
                                            if var_ic[c1[i]] > var_ic[c2[i]]: ic += var_ic[c2[i]]
                                            else: ic += var_ic[c1[i]]
                                        #X#print c1[i], var_ic[c1[i]], c2[i], var_ic[c2[i]]
                                    ## Setup Match Assessment counts of positions relationships ##
                                    assessment1 = {'Overlap':len(v1)-len(c1)}     # Dictionary of comtype and count for motif1
                                    assessment2 = {'Overlap':len(v2)-len(c2)}     # Dictionary of comtype and count for motif2
                                    for type in comtype.keys():
                                        assessment1[type] = 0
                                        for x in m1:
                                            if x in comtype[type]: assessment1[type] += 1
                                        assessment2[type] = 0
                                        for x in m2:
                                            if x in comtype[type]: assessment2[type] += 1
                                    #x#self.deBug('%s vs %s\n\t1: %s = %s\n\t2: %s = %s' % (c1,c2,m1,assessment1,m2,assessment2))
                                    ## Reject bad match ##
                                    if assessment1['Bad'] > self.stat['Mismatches'] or assessment2['Bad'] > self.stat['Mismatches']:    
                                        continue    # Incompatible match - bad positions
                                    if assessment1['Good'] < minshare or assessment2['Good'] < minshare:    
                                        continue    # Not enough matched positions of any kind
                                    if self.stat['MatchFix'] in [1,3] and assessment1['Fix'] < motif1.slimFix():  
                                        continue    # Not all fixed positions in motif1 matched
                                    if self.stat['MatchFix'] in [2,3] and assessment2['Fix'] < motif2.slimFix():  # Bad!
                                        continue    # Not all fixed positions in motif2 matched
                                    ## Update vmatch as appropriate ##
                                    vmatch = self._updateVMatch(vmatch,ic,assessment1,assessment2,(c1,c2))
                            ## Assess best vmatch and assign details to match ##
                            match = self._interpretMatch(vmatch)
                            
                    ## 4. Output match statistics ##
                    #x#self.deBug('%s v %s: %s' % (seq1,seq2,match))
                    if match:
                        outdict['Sim1'] = match[0]
                        outdict['Sim2'] = match[1]
                        outdict['MatchPos'] = match[2]
                        if motif1.stat['IC'] > motif2.stat['IC']: normic = match[3] / motif2.stat['IC']
                        else: normic = match[3] / motif1.stat['IC']
                        if normic < self.stat['NormCut']: continue     # Not a good enough match to retain
                        if match[0] == 'Duplicate': outdict['Match'] = outdict['Motif1']
                        else: outdict['Match'] = match[-1]        # match[4]  # Summary of match
                        outdict['NormIC'] = '%.3f' % normic
                        outdict['MatchIC'] = '%.3f' % match[3]
                        #x#outdict['Score'] = '%.3f' % (normic * match[3])    # NormIC x MatchIC
                        outdict['Score'] = '%.3f' % (normic * match[2])    # NormIC x MatchPos
                        mx += 1
                        rje.delimitedFileOutput(self,resfile,headers,delimit,outdict)
                        self.list['Matches'].append(rje.combineDict({},outdict))

            self.log.printLog('\r#COMP','Comparing %s & %s: %s comparisons made.' % (hlist[1][0],hlist[2][0],rje.integerString(compx)))
            self.log.printLog('\r#COMP','%s motif matches output to %s. Rank by "Score" column.' % (rje.integerString(mx),resfile))
            return
        except:
            self.log.errorLog('Comparimotif.motifCompare is banjaxed. Slap Rich till he fixes it.',quitchoice=True)
#########################################################################################################################
    def speedskip(self,seq1,seq2):  ### Speedily checks that seq1 and seq2 *may* be able to meet match requirements
        '''
        Speedily checks that seq1 and seq2 *may* be able to meet match requirements.
        << returns True if they should be skipped, else False.
        '''
        minshare = self.stat['MinShare']
        shared = 0
        for a in self.obj['SlimList'].list['Alphabet'][0:]:
            x = 0
            while seq1.count(a) > x and seq2.count(a) > x: x += 1
            shared += x
        if shared >= minshare: return False
        return True
#########################################################################################################################
    def _positionRelationship(self,c1,c2):  ### Compares motif elements c1 and c2 and returns their relationship
        '''
        Compares motif elements c1 and c2 and returns their relationship. Each should be a fixed position, an ambiguous
        position, or a wildcard 'X'.
        >> c1:str = Element of variant of Motif1 in motif-motif comparison
        >> c2:str = Element of variant of Motif2 in motif-motif comparison
        << rel:tuple = relationship of (c1 vs c2, c2 vs c1)
        '''
        try:
            ### Obvious relationships first ###
            if c1 in ['.','X'] and c2 in ['.','X']: return ('wild','wild')      # Both wildcards
            elif c1 in ['$','^'] and c2 != c1: return ('mismatch','mismatch')
            elif c2 in ['$','^'] and c2 != c1: return ('mismatch','mismatch')     
            elif c1 in ['.','X']: return ('wilddegen','wildvar')          # One wildcard (matches anything)
            elif c2 in ['.','X']: return ('wildvar','wilddegen')          # One wildcard (matches anything)
            elif c1 == c2 and len(c1) == 1: return ('match','match')# Non-wildcard match
            elif c1 == c2 or rje.strSort(c1) == rje.strSort(c2): return ('ambmatch','ambmatch') # Ambiguity match
            elif len(c2) < 2 and c1.find(c2) >= 0: return ('degen','variant')       # Degenerate version
            elif len(c1) < 2 and c2.find(c1) >= 0: return ('variant','degen')       # Degenerate version
            elif self.ambVariant(c1,c2): return ('ambdegen','ambvar')       # Ambiguity degeneracy
            elif self.ambVariant(c2,c1): return ('ambvar','ambdegen')       # Ambiguity degeneracy
            ### Complex relationships ###                
            elif len(c1) > 1 and len(c2) > 1:                       # Two ambiguities
                ax = False
                for a in c1:
                    if c2.find(a) >= 0: return ('overlap','overlap')# Ambiguities overlap
                return ('badamb','badamb')                          # Ambiguities don't overlap
            elif len(c1) > 1: return ('badamb','ambmismatch')       # Ambiguity mismatch
            elif len(c2) > 1: return ('ambmismatch','badamb')       # Ambiguity mismatch
            ### No relationship ###
            else: return ('mismatch','mismatch')                    # Fixed mismatch
        except:
            self.log.errorLog('Run-time error in CompariMotif._positionRelationship()')
            raise
#########################################################################################################################
    def ambVariant(self,el1,el2):   ### Returns whether el2 is an ambiguous variant of el1
        '''Returns whether el2 is an ambiguous variant of el1.'''
        for a in el2:
            if not el1.find(a) >= 0: return False
        return True     # All AAs in el2 are found in el1
#########################################################################################################################
    def _updateVMatch(self,vmatch,ic,assessment1,assessment2,overlap):  ### Assesses and updates current best variant match stats
        '''
        Assesses and updates current best variant match stats.
        >> vmatch:tuple = current best match: (ic,assessment1,assessment2)
        >> ic:float = Information Content of current match
        >> assessment1:dict = counts of different match types for motif1
        >> assessment2:dict = counts of different match types for motif2
        << newmatch:tuple = updated best match: (ic,assessment1,assessment2)
        '''
        try:
            ### Check whether existing vmatch is best###
            if vmatch and ic < vmatch[0]: return vmatch
            elif vmatch and ic == vmatch[0]:
                if vmatch[1]['Good'] > assessment1['Good']: return vmatch
                elif vmatch[1]['Good'] == assessment1['Good'] and vmatch[1]['Fix'] > assessment1['Fix']: return vmatch
                elif vmatch[1]['Good'] == assessment1['Good'] and vmatch[1]['Fix'] == assessment1['Fix'] and ic <= vmatch[0]: return vmatch
            ### Return new vmatch ###
            newmatch = (ic,{},{},overlap)
            for key in assessment1.keys():
                newmatch[1][key] = assessment1[key]
                newmatch[2][key] = assessment2[key]
            return newmatch
        except:
            self.log.errorLog('Run-time error in CompariMotif._updateVMatch()')
            raise
#########################################################################################################################
    def _interpretMatch(self,vmatch):   ### Translates vmatch stats into match information for output
        '''
        Translates vmatch stats into match information for output.
        >> vmatch:tuple = current best match: (ic,assessment1,assessment2,overlap tuple)
        << match:list = [matchtype1,matchtype2,positions matched,information content]
        '''
        try:
            ### No match if no vmatch ###
            if not vmatch: return []
            ### Setup basic match stats ###
            match = ['','',vmatch[1]['Good'],vmatch[0],vmatch[3]]
            #X#self.deBug(match)
            ### Define match types ###
            for m in [1,2]:
                if vmatch[m]['Degen'] == 0 and vmatch[m]['Variant'] == 0: match[m-1] = 'Exact'
                elif vmatch[m]['Degen'] > 0 and vmatch[m]['Variant'] > 0: match[m-1] = 'Complex'
                elif vmatch[m]['Degen'] > 0: match[m-1] = 'Degenerate'
                elif vmatch[m]['Variant'] > 0: match[m-1] = 'Variant'
                else: match[m-1] = 'Problem'
            ### Define match lengths
                if vmatch[1]['Overlap'] == 0 and vmatch[2]['Overlap'] == 0: match[m-1] += ' Match'
                elif vmatch[1]['Overlap'] > 0 and vmatch[2]['Overlap'] > 0: match[m-1] += ' Overlap'
                elif vmatch[m]['Overlap'] == 0: match[m-1] += ' Subsequence'
                elif vmatch[m]['Overlap'] > 0: match[m-1] += ' Parent'
                else: match[m-1] += ' Problem'
            ### Generate description of overlapping portion ###
            (c1,c2) = vmatch[3]
            common = []
            for i in range(len(c1)):
                if c1[i] == c2[i]: common.append(c1[i])
                elif c1[i] == '.': common.append(c2[i].lower())
                elif c2[i] == '.': common.append(c1[i].lower())
                else:
                    common.append('')
                    for aa in c1[i].lower():
                        if c2[i].lower().find(aa) >= 0: common[-1] += aa
                    if not common[-1]: common[-1] = '*'
                clow = common[-1] == common[-1].lower()
                if len(common[-1]) > 15:
                    common[-1] = '^%s' % string.join(list(sets.Set(rje_slim.default_aas).difference(common[-1].upper())),'')
                    if clow: common[-1] = common[-1].lower()
                if len(common[-1]) > 1: common[-1] = '[%s]' % common[-1]
            match[4] = string.join(common,'')
            ### Return ###
            return match
        except:
            self.log.errorLog('Run-time error in CompariMotif._interpretMatch()')
            raise
#########################################################################################################################
    ### <4> ### Additional outputs                                                                                      #
#########################################################################################################################
    def compareXGMML(self,save=True):     ### Outputs XGMML format file for Cytoscape
        '''
        Outputs XGMML format file for Cytoscape.
        >> save:bool = whether to save file (True) or return XGMML object (False) [True]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['XGMML']: return None
            xgmml = rje_xgmml.XGMML(self.log,self.cmd_list)
            motifs = self.obj['SlimList']
            compdb = self.obj['SearchDB']
            motiflist = motifs.motifs() + compdb.motifs()
            if motifs.info['SearchDB'] == motifs.info['Motifs']: motiflist = motifs.motifs()[0:]
            ## ~ [1a] Node attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodeatt = xgmml.dict['NodeAtt'] = {'Motif':'string','Pattern':'string','Description':'string',
                                               'PosLength':'real','MinLength':'real','MaxLength':'real',
                                               'FixLength':'real','IC':'real','Dataset':'string'}
            for cyt in ['shape','fillColor','size','fontSize','borderColor','font']: xgmml.dict['NodeAtt']['node.%s' % cyt] = 'string' 
            edgeatt = xgmml.dict['EdgeAtt'] = {'NormIC':'real','Score':'real','MatchPos':'real','MatchIC':'real','Match':'string'}
            for cyt in ['lineStyle','color','sourceArrowShape','sourceArrowColor','targetArrowShape','targetArrowColor','lineWidth']: xgmml.dict['EdgeAtt']['edge.%s' % cyt] = 'string' 
                
            ### ~ [2] Output XGMML files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xgmml.info['Name'] = rje.baseFile(self.info['Resfile'])
            xfile = '%s.xgmml' % xgmml.info['Name']

            ## ~ [2b] Update node data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodes = xgmml.dict['Node'] = {}
            for motif in motiflist:
                name = motif.info['Name']
                nodes[name] = {'Motif':motif.info['Name'],'Pattern':motif.pattern(),
                               'Description':motif.info['Description'],'IC':'%.3f' % motif.stat['IC'],
                               'MaxLength':motif.slimLen(),'MinLength':motif.slimMinLen(),'FixLength':motif.slimFix(),
                               'PosLength':motif.slimPos()}
                if motif in compdb.motifs():
                    nodes[name]['node.shape'] = 'rectangle'
                    nodes[name]['node.fillColor'] = '255,153,153'
                    nodes[name]['Dataset'] = compdb.info['Name']
                else:
                    nodes[name]['node.shape'] = 'ellipse'
                    nodes[name]['node.fillColor'] = '153,153,255'
                    nodes[name]['Dataset'] = motifs.info['Name']
                nodes[name]['node.font'] = 'Impact,plain,10'
                nodes[name]['node.fontSize'] = '10.0'
                nsize = 25 + (motif.stat['IC'] *5)
                nodes[name]['node.size'] = '%d.0' % nsize

            ## ~ [2c] Update edge data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xmatch = []
            edges = xgmml.dict['Edge'] = {}
            for match in self.list['Matches']:
                m1 = match['Name1']
                m2 = match['Name2']
                if match['Sim1'] == 'Duplicate': s1 = 'dup'
                else:
                    s1 = string.split(match['Sim1'].lower())
                    s1 = s1[0][0] + s1[1][0]
                if match['Sim2'] == 'Duplicate': s2 = 'dup'
                else:
                    s2 = string.split(match['Sim2'].lower())
                    s2 = s2[0][0] + s2[1][0]
                if s1 not in edges: edges[s1] = {}
                if (m2,m1) in xmatch: continue   # Reverse match already there
                xmatch.append((m1,m2))
                edges[s1][(m1,m2)] = {}
                for result in ['NormIC','Score','MatchPos','Match','MatchIC']: edges[s1][(m1,m2)][result] = match[result]
                edges[s1][(m1,m2)]['edge.lineStyle'] = 'SOLID'
                if s1[0] == 'c': edges[s1][(m1,m2)]['edge.lineStyle'] = 'LONG_DASH'
                edges[s1][(m1,m2)]['edge.color'] = '0,0,0'
                edges[s1][(m1,m2)]['edge.lineWidth'] = '%.2f' % (1 + string.atof(match['NormIC']))
                if s1[1] in ['s','o']:
                    edges[s1][(m1,m2)]['edge.targetArrowShape'] = 'DELTA'
                    edges[s1][(m1,m2)]['edge.targetArrowColor'] = edges[s1][(m1,m2)]['edge.color']
                if s1[1] in ['p','o']:
                    edges[s1][(m1,m2)]['edge.sourceArrowShape'] = 'DELTA'
                    edges[s1][(m1,m2)]['edge.sourceArrowColor'] = edges[s1][(m1,m2)]['edge.color']
            
            ## ~ [2d] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if save: xgmml.saveXGMML(xfile)
            return xgmml
        except: self.log.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
### End of SECTION II: CompariMotif Class                                                                               #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try: CompariMotif(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.printLog('#END','User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
