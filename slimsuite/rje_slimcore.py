#!/usr/bin/python

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
Module:       rje_slimcore
Description:  Core module/object for SLiMFinder and SLiMSearch
Version:      1.7
Last Edit:    03/06/10
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is primarily to contain core dataset processing methods for both SLiMFinder and SLiMSearch to inherit and
    use. This primarily consists of the options and methods for masking datasets and generating UPC. This module can
    therefore be run in standalone mode to generate UPC files for SLiMFinder or SLiMSearch.

    In addition, the secondary MotifSeq and Randomise functions are handled here.
    
Secondary Functions:
    The "MotifSeq" option will output fasta files for a list of X:Y, where X is a motif pattern and Y is the output file.

    The "Randomise" function will take a set of input datasets (as in Batch Mode) and regenerate a set of new datasets
    by shuffling the UPC among datasets. Note that, at this stage, this is quite crude and may result in the final
    datasets having fewer UPC due to common sequences and/or relationships between UPC clusters in different datasets.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options ### 
    seqin=FILE      : Sequence file to search [None]
    batch=LIST      : List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    maxseq=X        : Maximum number of sequences to process [500]
    maxupc=X        : Maximum UPC size of dataset to process [0]
    sizesort=X      : Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
    walltime=X      : Time in hours before program will abort search and exit [1.0]
    resdir=PATH     : Redirect individual output files to specified directory (and look for intermediates) [SLiMFinder/]
    buildpath=PATH  : Alternative path to look for existing intermediate files [SLiMFinder/]
    force=T/F       : Force re-running of BLAST, UPC generation and SLiMBuild [False]
    dna=T/F         : Whether the sequences files are DNA rather than protein [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Evolutionary Filtering Options ###
    efilter=T/F     : Whether to use evolutionary filter [True]
    blastf=T/F      : Use BLAST Complexity filter when determining relationships [True]
    blaste=X        : BLAST e-value threshold for determining relationships [1e=4]
    altdis=FILE     : Alternative all by all distance matrix for relationships [None]
    gablamdis=FILE  : Alternative GABLAM results file [None] (!!!Experimental feature!!!)
    domtable=FILE   : Domain table containing domain ("Type") and sequence ("Name") pairings for additional UPC [None]
    homcut=X        : Max number of homologues to allow (to reduce large multi-domain families) [0]
    extras=T/F      : Whether to generate additional output files (distance matrices etc.) [True]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Input Masking and AA Frequency Options ###
    masking=T/F     : Master control switch to turn off all masking if False [True]
    dismask=T/F     : Whether to mask ordered regions (see rje_disorder for options) [False]
    consmask=T/F    : Whether to use relative conservation masking [False]
    ftmask=LIST     : UniProt features to mask out [EM,DOMAIN,TRANSMEM]
    imask=LIST      : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
    fakemask=T/F    : Whether to invoke UniFake to generate additional features for masking [False]
    compmask=X,Y    : Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    casemask=X      : Mask Upper or Lower case [None]
    metmask=T/F     : Masks the N-terminal M (can be useful if SLiMFinder termini=T) *Also maskm=T/F* [True]
    posmask=LIST    : Masks list of position-specific aas, where list = pos1:aas,pos2:aas *Also maskpos=LIST* [2:A]
    aamask=LIST     : Masks list of AAs from all sequences (reduces alphabet) []
    motifmask=X     : List (or file) of motifs to mask from input sequences []
    logmask=T/F     : Whether to output the log messages for masking of individual sequences to screen [False]
    masktext=X      : Text ID to over-ride automated masking text and identify specific masking settings [None]
    maskpickle=T/F  : Whether to save/load pickle of masked input data, independent of main pickling [False]
    maskfreq=T/F    : Whether to use masked AA Frequencies (True), or (False) mask after frequency calculations [True]
    aafreq=FILE     : Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]    
    smearfreq=T/F   : Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    qregion=X,Y     : Mask all but the region of the query from (and including) residue X to residue Y [0,-1]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Advanced Output Options ###
    targz=T/F       : Whether to tar and zip dataset result files (UNIX only) [False]
    pickle=T/F      : Whether to save/use pickles [True]
    savespace=0     : Delete "unneccessary" files following run (best used with targz): [0]
                        - 0 = Delete no files
                        - 1 = Delete all bar *.upc, *.pickle and *.occ.csv files (Pickle excluded from tar.gz with this setting)
                        - 2 = Delete all bar *.upc, *.pickle files (Pickle excluded from tar.gz with this setting)
                        - 3 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Additional Functions I: MotifSeq ###    
    motifseq=LIST   : Outputs fasta files for a list of X:Y, where X is the pattern and Y is the output file []
    slimbuild=T/F   : Whether to build motifs with SLiMBuild. (For combination with motifseq only.) [True]

    ### Additional Functions II: Randomised datasets ###
    randomise=T/F   : Randomise UPC within batch files and output new datasets [False]
    randir=PATH     : Output path for creation of randomised datasets [Random/]
    randbase=X      : Base for random dataset name [rand]
    randsource=FILE : Source for new sequences for random datasets (replaces UPCs) [None]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Uses general modules: copy, glob, math, os, string, sys, time
Uses RJE modules: rje, rje_blast, rje_dismatrix_V2, rje_seq, rje_scoring, rje_slim, rje_slimlist
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, math, os, pickle, random, string, sys, time
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_blast, rje_motiflist, rje_seq, rje_sequence, rje_scoring, rje_slim, rje_slimcalc, rje_slimlist
import rje_motif_V3 as rje_motif
import rje_dismatrix_V2 as rje_dismatrix
import comparimotif_V3 as comparimotif
import unifake
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on SLiMFinder 3.1.
    # 0.1 - Tidied with respect to SLiMFinder and SLiMSearch.
    # 0.2 - Added DNA mode.
    # 0.3 - Added relative conservation masking.
    # 0.4 - Altered TarGZ to *not* include *.pickle.gz
    # 1.0 - Standardised masking options. Added motifmask and motifcull.
    # 1.1 - Checked/updated Randomise option.
    # 1.2 - Added generation of UniFake file to maskInput() method. (fakemask=T)
    # 1.3 - Added aamask effort
    # 1.4 - Added DomTable option
    # 1.5 - Added masktext and maskpickle options to accelerate runs with large masked datasets.
    # 1.6 - Fixed occurrence table bugs.
    # 1.7 - Added SizeSort and NewUPC. Add server #END statements.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Reduce methods and attributes of SLiMCore object to only those necessary for function.
    # [ ] : Reduce commandline options accordingly - proper use of coreCmd() and add coreDefaults().
    # [ ] : Add proper use of seqbase and basefile for SLiMSearch results output.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('SLiMCore', '1.7', 'April 2010', '2007')
    description = 'SLiMFinder/SLiMSearch Core Module'
    author = 'Richard J. Edwards'
    comments = ['Please report bugs to R.Edwards@Soton.ac.uk']
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
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)
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
### CONSTANTS ###                                                                                                     
wildcards = ['.','X','x']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMCore Class                                                                                          #
#########################################################################################################################
class SLiMCore(rje.RJE_Object):     
    '''
    SLiMFinder Class. Author: Rich Edwards (2007).

    Info:str
    - AAFreq = Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    - AltDis = Alternative all by all distance matrix for relationships [None]
    - Build = String given summary of key SLiMBuild options
    - BuildPath = Alternative path to look for existing intermediate files [SLiMFinder/]
    - CaseMask = Mask Upper or Lower case [None]
    - CompMask = Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    - DomTable = Domain table containing domain ("Type") and sequence ("Name") pairings for additional UPC [None]
    - GablamDis = Alternative GABLAM results file [None]
    - Input = Original name (and path) of input file
    - MaskText = Text ID to over-ride automated masking text and identify specific masking settings [None]
    - MotifMask = List (or file) of motifs to mask from input sequences []
    - RanDir = Output path for creation of randomised datasets [./]
    - Randbase = Base for random dataset name [rand_]
    - RandSource = Source for new sequences for random datasets (replaces UPCs) [None]
    - ResDir = Redirect individual output files to specified directory [SLiMFinder/]
    
    Opt:boolean
    - ConsMask = Whether to use relative conservation masking [False]
    - DisMask = Whether to mask ordered regions (see rje_disorder for options) [False]
    - DNA = Whether the sequences files are DNA rather than protein [False]
    - Force = whether to force recreation of key files [False]
    - EFilter = Whether to use evolutionary filter [True]
    - Extras = Whether to generate additional output files (alignments etc.) [False]
    - FakeMask = Whether to invoke UniFake to generate additional features for masking [False]
    - LogMask = Whether to log the masking of individual sequences [True]
    - Masked = Whether dataset has been masked [False]
    - Masking = Master control switch to turn off all masking if False [True]
    - MaskM = Masks the N-terminal M (can be useful if termini=T) [False]
    - MaskPickle = Whether to save/load pickle of masked input data, independent of main pickling [False]
    - OccStatsCalculated = Whether OccStats have been calculated for all occurrence [False]
    - MaskFreq = Whether to mask input before any analysis, or after frequency calculations [True]
    - Pickle = Whether to save/use pickles [True]
    - Randomise = Randomise UPC within batch files [False]
    - SlimBuild = Whether to build motifs with SLiMBuild. (For combination with motifseq only.) [True]
    - SmearFreq = Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    - TarGZ = Whether to tar and zip dataset result files (UNIX only) [False]
    - Test = Special Test parameter for experimentation with code [False]
    - Webserver = Generate additional webserver-specific output [False]

    Stat:numeric
    - HomCut = Max number of homologues to allow (to reduce large multi-domain families) [0]
    - MaxSeq = Maximum number of sequences to process [500]
    - MaxUPC = Maximum UPC size of dataset to process [0]
    - MST = MST corrected size for whole dataset
    - SaveSpace = Delete "unneccessary" files following run (see Manual for details) [0]
    - SizeSort = Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
    - StartTime = Starting time in seconds (for output when using shared log file)
    - WallTime = Time in hours before program will terminate [1.0]

    List:list
    - AAMask = Masks list of AAs from all sequences (reduces alphabet) []
    - Alphabet = List of characters to include in search (e.g. AAs or NTs) 
    - Batch = List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    - FTMask = UniProt features to mask out [EM,DOMAIN,TRANSMEM]
    - Headers = Headers for main SLiMFinder output table
    - IMask = UniProt features to inversely ("inclusively") mask [IM]
    - QRegion = Mask all but the region of the query from (and including) residue X to residue Y [0,-1]
    - SigSlim = List of significant SLiMs - matches keys to self.dict['Slim(Freq)'] - *in rank order*
    - UP = List of UP cluster tuples
    - Warning = List of text (log) warnings to reproduce at end of run
     
    Dict:dictionary
    - AAFreq = AA frequency dictionary for each seq / UPC
    - DimFreq = Frequency of dimers of each X length per upc {upc:[freqs]}
    - Dimers = main nested dictionary for SLiMFinder {Ai:{X:{Aj:{'UP':[UPC],'Occ':[(Seq,Pos)]}}}}
    - ElementIC = dictionary of {Motif Element:Information Content}
    - Extremf. = Dictionary of {length:extremferroni correction}
    - MaskPos = Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    - MotifSeq = Dictionary of {pattern:output file for sequences}
    - MST = MST corrected size for UPC {UPC:MST}
    - Slim = main dictionary containing SLiMs with enough support {Slim:{'UPC':[UPC],'Occ':[Seq,Pos]}}
    - SeqOcc = dictionary of {Slim:{Seq:Count}}
    - SmearFreq = Smeared AA frequency dictionary for IC calculations

    Obj:RJE_Objects
    - SeqList = main SeqList object containing dataset to be searched
    - SlimList = MotifList object handling motif stats and filtering options
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AAFreq','Input','CaseMask','CompMask','ResDir','RanDir','Randbase','Build','BuildPath',
                         'AltDis','GablamDis','RandSource','DomTable','SeqIn','MaskText']
        self.optlist = ['BuildProb','DisMask','Force','MaskFreq','Test','LogMask','MaskM','Masked','Masking','EFilter',
                        'SlimBuild','Randomise','Extras','TarGZ','SmearFreq','OccStatsCalculated','DNA','ConsMask',
                        'FakeMask','MaskPickle','Webserver']
        self.statlist = ['StartTime','MST','WallTime','MaxSeq','SaveSpace','MaxUPC','HomCut','Extras','SizeSort']
        self.listlist = ['AAMask','Alphabet','Batch','FTMask','IMask','SigSlim','UP','Headers','Warning','QRegion']
        self.dictlist = ['AAFreq','DimFreq','Dimers','Slim','SeqOcc','Extremf.','ElementIC','MST','MotifSeq','MaskPos',
                         'SmearFreq']
        self.objlist = ['SeqList','SlimList']
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.coreDefaults()
#########################################################################################################################
    def coreDefaults(self): ### Sets the default parameters used by SLiMCore functions (for easy inheritance).
        '''### Sets the default parameters used by SLiMCore functions (for easy inheritance).'''
        self.setInfo({'CompMask':'5,8','ResDir':rje.makePath('SLiMCore/'),'Input':'','Basefile':'','MotifMask':'',
                      'RanDir':rje.makePath('Random/'),'Randbase':'rand','BuildPath':rje.makePath('SLiMCore/'),
                     'MaskText':'None'})
        self.setStat({'StartTime':time.time(),'WallTime':1.0,'MaxSeq':500,'SaveSpace':0,'MaxUPC':0,'HomCut':0,'Extras':1,
                      'SizeSort':0})
        self.setOpt({'MaskFreq':True,'Extras':True,'Masked':False,'MaskM':True,'Randomise':False,'TarGZ':False,
                     'Force':False,'SmearFreq':False,'DisMask':False,'Test':False,'OccStatsCalculated':False,
                     'LogMask':False,'DNA':False,'ConsMask':False,'Pickle':True,'FakeMask':False,'MaskPickle':False,
                     'Webserver':False})
        self.list['FTMask'] = ['EM','Domain','TRANSMEM']
        self.list['Batch'] = ['*.dat','*.fas']
        self.dict['MaskPos'] = {2:'A'}
        self.dict['SmearFreq'] = {}
        ### Other Attributes ###
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
        self.obj['SlimList'].obj['SLiMCalc'].setupFilters(slimheaders=[],occheaders=[])    ### Sets up SLiM/Occ Filters
#########################################################################################################################
    def _cmdList(self): self.coreCmd()
#########################################################################################################################
    def coreCmd(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'file',['AAFreq','AltDis','GablamDis','RandSource','DomTable','SeqIn'])
                self._cmdReadList(cmd,'info',['CaseMask','CompMask','Randbase','MotifMask','MaskText'])
                self._cmdReadList(cmd,'path',['ResDir','RanDir','BuildPath'])
                self._cmdReadList(cmd,'int',['SaveSpace','MaxSeq','MaxUPC','HomCut','Extras','SizeSort'])
                self._cmdReadList(cmd,'stat',['WallTime'])
                self._cmdReadList(cmd,'opt',['DisMask','Force','MaskFreq','Extras','Test','LogMask','Randomise','MaskM',
                                             'EFilter','TarGZ','Masking','SlimBuild','SmearFreq','DNA','ConsMask',
                                             'Pickle','FakeMask','MaskPickle','Webserver'])
                self._cmdRead(cmd,'opt','MaskM','metmask')
                self._cmdReadList(cmd,'list',['FTMask','IMask','AAMask','Alphabet','Batch','QRegion'])
                self._cmdReadList(cmd,'cdict',['MotifSeq','MaskPos'])
                self._cmdRead(cmd,'cdict','MaskPos','posmask')
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self.setAlphabet()
        self.setQRegion()
#########################################################################################################################
    ### <2> ### Simple Stats Methods                                                                                    #
#########################################################################################################################
    def occFilter(self): return self.obj['SlimList'].obj['SLiMCalc'].dict['OccFilter']
    def statFilter(self): return self.obj['SlimList'].obj['SLiMCalc'].dict['SlimFilter']
    def occStats(self): return self.obj['SlimList'].obj['SLiMCalc'].list['SLiMCalc']
    def dataset(self): return os.path.split(self.info['Basefile'])[1]
    def seqs(self): return self.obj['SeqList'].seq[0:]
    def seqNum(self): return self.obj['SeqList'].seqNum()
    def units(self): return self.obj['SeqList'].units()
    def aaNum(self): return self.obj['SeqList'].aaTotal()
    def nonX(self): return self.obj['SeqList'].aaTotal(nonx=True)
#########################################################################################################################
    def slimNum(self): return len(self.dict['Slim'])
#########################################################################################################################
    def slimUP(self,slim):  ### Returns UP Num for SLiM if in dictionary, else 0.
        if self.dict['Slim'].has_key(slim): return len(self.dict['Slim'][slim]['UP'])
        return 0
#########################################################################################################################
    def slimOccNum(self,slim,upc=None):    ### Returns number of occ of Slim in given UPC
        '''Returns number of occ of Slim in given UPC.'''
        ## No UPC - return total! ##
        if not upc: return len(self.dict['Slim'][slim]['Occ'])
        ## See if SLiM occurs at all! ##
        if not self.dict['SeqOcc'].has_key(slim): return 0
        ## Mean over UPC ##
        sx = 0
        for seq in upc:
            if self.dict['SeqOcc'][slim].has_key(seq) and self.dict['SeqOcc'][slim][seq] > sx:
                sx = self.dict['SeqOcc'][slim][seq]
        return sx
#########################################################################################################################
    def slimSeqNum(self,slim):  ### Returns the number of sequences SLiM occurs in.
        '''Returns the number of sequences SLiM occurs in.'''
        if not self.dict['SeqOcc'].has_key(slim): return 0      ## See if SLiM occurs at all! ##
        return len(self.dict['SeqOcc'][slim])
#########################################################################################################################
    def UPNum(self): return len(self.list['UP'])
#########################################################################################################################
    def getUP(self,seq):
        for upc in self.list['UP']:
            if seq in upc: return upc
        return None
#########################################################################################################################
    def slimIC(self,slim,usefreq=False):   ### Returns IC of slim  #!# Add aafreq & wildcard pen etc? 
        '''Returns IC of slim. Does not account for variable length wildcards!'''
        ic = 0.0
        slist = string.split(slim,'-')
        if usefreq and ('SmearFreq' not in self.dict or not self.dict['SmearFreq']): self.dict['SmearFreq'] = self.smearAAFreq(False)
        i = 0
        while i < (len(slist)):     #!# Update this for new calculation and DNA?! #!#
            el = slist[i]
            if not self.dict['ElementIC'].has_key(el):
                #x#self.dict['ElementIC'][el] = rje_motif.elementIC(el,callobj=self)
                if usefreq: self.dict['ElementIC'][el] = rje_slim.weightedElementIC(el,self.dict['SmearFreq'])
                elif self.opt['DNA']: self.dict['ElementIC'][el] = rje_slim.weightedElementIC(el,rje_slim.even_ntfreq)
                else: self.dict['ElementIC'][el] = rje_slim.weightedElementIC(el,rje_slim.even_aafreq)
            ic += self.dict['ElementIC'][el]
            i += 2
        return ic
#########################################################################################################################
    ### <3> ### General Run Methods                                                                                     #
#########################################################################################################################
    def run(self,batch=False):  ### Main Run Method. SLiMFinder, SLiMSearch etc. should have their own methods.
        '''
        Main Run Method. SLiMFinder, SLiMSearch etc. should have their own methods.
        1. Check for randomise function and execute if appropriate
        2. Input:
            - Read sequences into SeqList
            - or - Identify appropriate Batch datasets and rerun each with batch=True
        3. Masking and UPC construction.
        4. MotifSeq option if desired.
        >> batch:bool [False] = whether this run is already an individual batch mode run.
        '''
        try:### ~ [1] Randomise Function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Randomise'] and not batch: return self.randomise()

            ### ~ [2] Input/Batch handling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            self.setupBasefile()
            ## ~ [2a] Batch Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not batch and self.seqNum() < 1:   # No sequences loaded - use batch mode
                #batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
                #self.log.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
                batchfiles = self.batchFiles()
                if not batchfiles: self.log.errorLog('No input files found!',printerror=False)
                else:
                    mycmd = self.cmd_list[0:]
                    self.list['Batch'] = []
                    bx = 0
                    for infile in batchfiles:
                        bx += 1
                        self.log.printLog('#BATCH','Batch running %s' % infile)
                        bsf = self.newBatchRun(infile)
                        bsf.run(batch=True)
                        self.log.printLog('#BATCH','Batch file %s run. Cleaning up for next file.' % infile)
                        del bsf.obj
                        del bsf.list
                        del bsf.dict
                        del bsf
                        self.opt['Append'] = True
                        self.log.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(bx),rje.integerString(len(batchfiles)-bx)),log=False)
                if self.opt['Win32'] and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return
            ## ~ [2b] Check whether to bother running dataset at all - Check Input versus Min and Max Seq ~ ##
            if self.stat['MaxSeq'] > 0 and self.stat['MaxSeq'] < self.obj['SeqList'].seqNum():
                self.log.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.integerString(self.obj['SeqList'].seqNum()),rje.integerString(self.stat['MaxSeq'])))
                return False

            ### ~ [3] UPC and Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.stat['StartTime'] = time.time()
            ## ~ [3a] UPC and MinOcc settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.makeUPC(): return self.log.errorLog('Error during %s.makeUPC(). Abandoning %s run' % (self.prog(),self.dataset()),printerror=False)
            if self.stat['MaxUPC'] > 0 and self.stat['MaxUPC'] < self.UPNum():
                self.log.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). Run aborted.' % (self.UPNum(),self.stat['MaxUPC']))
                return
            ## ~ [3b] AA Frequency Calculations made early as needed superficially in SLiMBuild ~~~ ##
            if self.dict['MotifSeq']:
                self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
                if self.opt['Masked']: self.obj['SeqList'].saveFasta(seqfile='%s.%s.masked.fas' % (self.info['Basefile'],self.maskText()))

            ### ~ [4] Special MotifSeq Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                self.motifSeq()
    
            ### End ###
                self.printLog('#CORE','SLiMCore results output to %s*' % (self.info['Basefile']))
            if self.opt['Win32'] and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            return True
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.stat['WallTime'] <= 0 or (time.time() - self.stat['StartTime']) < (self.stat['WallTime']*3600): raise
            if self.list['Headers']: self.results(aborted=True)
            return False # Walltime reached
        except:
            self.errorLog('Error in %s.run()' % self.prog(),printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def batchFiles(self,pickup=[]):   ### Returns batch file list
        '''Returns batch file list.'''
        ### ~ [0] ~ Standard Batch file list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
        self.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
        try:
            if not self.stat['SizeSort']: return batchfiles
        except: return batchfiles
        ### ~ [1] SizeSort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        sizedict = {}
        for file in batchfiles:
            if os.path.split(rje.baseFile(file))[1] in pickup:
                self.printLog('#PICKUP','Skipping %s batch file %s' % (self.info['RunID'],file))
                continue
            fsize = os.path.getsize(file)
            if fsize not in sizedict: sizedict[fsize] = []
            sizedict[fsize].append(file)
        sizes = rje.sortKeys(sizedict,revsort=self.stat['SizeSort']<0)
        newbatch = []       # New list of batch files in size order, starting with the biggest
        for fsize in sizes: newbatch += sizedict.pop(fsize)
        return newbatch
#########################################################################################################################
    def wallTime(self):     ### Exits if walltime has been reached
        '''Exits if walltime has been reached.'''
        if self.stat['WallTime'] <= 0 or (time.time() - self.stat['StartTime']) < (self.stat['WallTime']*3600): return
        self.errorLog('%s Walltime (%.2f hours) reached! Try increasing WallTime=X.' % (self.prog(),self.stat['WallTime']),printerror=False)
        self.serverEnd('Walltime',exit=False)
        sys.exit()
#########################################################################################################################
    def newBatchRun(self,infile):   ### Returns SLiMFinder object for new batch run
        '''Returns SLiMCore object for new batch run.'''
        return SLiMCore(self.log,self.cmd_list[0:] + ['seqin=%s' % infile,'append=%s' % self.opt['Append']])
#########################################################################################################################
    def serverEnd(self,endcause,details=None,exit=True):  ### Adds #END statement for webserver
        '''
        Adds #END statement for webserver.
        >> endcause:str = Identifier for ending program run.
        >> details:str [None] = Extra information for certain end causes only.
        >> exit:bool [True] = whether to terminate program after #END statement.
        '''
        ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:
            dtxt = 'Download install a local copy of SLiMFinder <a href="http://www.southampton.ac.uk/~re1u06/software/slimsuite/">here</a> for larger runs.'
            if not self.opt['Webserver']:
                if endcause == 'Crash': sys.exit(1)
                else: return
            ### ~ [2] END statements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if endcause == 'Walltime': end = '%s server Walltime (%.2f hours) reached! Suggestions: reduce dataset size, increase masking, reduce slimlen and/or maxwild, switch ambiguity off.' % (self.prog(),self.stat['WallTime'])
            elif endcause == 'MaxSeq': end = 'Server limit on unrelated sequence (maxseq=%s) exceeded: %s input seq. Analysis terminated.' % (rje.integerString(self.stat['MaxSeq']),rje.integerString(self.seqNum()))
            elif endcause == 'MaxSize': end = 'Server limit on size of dataset (maxsize=%s %s) exceeded: %s input seq. Analysis terminated.' % (rje.integerString(self.stat['MaxSize']),self.units(),rje.integerString(self.aaNum()))
            elif endcause == 'FewSeq': end = 'Insufficient sequences for minimum SLiM occurrence setting (minocc=%d): %d input sequences; %d+ unrelated sequences needed. Run aborted.' % (self.stat['MinOcc'],self.seqNum(),self.stat['MinOcc'])
            elif endcause == 'MaxUPC': end = 'Server limit on unrelated sequence (maxupc=%s) exceeded: %s unrelated input sequences (UPC). Analysis terminated.' % (rje.integerString(self.stat['MaxUPC']),rje.integerString(self.UPNum()))
            elif endcause == 'FewUPC': end = 'Insufficient <i>unrelated</i> sequences for minimum SLiM occurrence setting (minocc=%d): %d unrelated input sequences (UPC); %d+ unrelated sequences needed. Run aborted.' % (self.stat['MinOcc'],self.UPNum(),self.stat['MinOcc'])
            elif endcause == 'NoSLiM': end = 'No SLiMs for output. Try relaxing probability threshold (probcut) and setting number of SLiMs to return (topranks) > 0.'
            elif endcause == 'NoMotif': end = 'No SLiMs retained for searching. Check format of input file and information content of motifs.'
            elif endcause == 'Crash':
                if details: end = 'Sorry - %s server crashed during %s! Please contact webmaster with job ID.' % (self.prog(),details)
                else: end = 'Sorry - %s server crashed! Please contact webmaster with job ID.' % self.prog()
            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else: end = '%s server run ended. See log file for details.' % self.prog()
            if endcause in ['Walltime','MaxSeq','MaxUPC']: end = '%s %s' % (end,dtxt)
            self.printLog('#END',end)
            if exit:
                if endcause == 'Crash': sys.exit(1)
                else: sys.exit(0)
        except SystemExit: raise
        except: self.errorLog('ServerEnd error.')
#########################################################################################################################
    ### <4> ### Setup/Input Methods                                                                                     #
#########################################################################################################################
    def setupBasefile(self):    ### Sets up self.info['Basefile'] and self.info['Input']
        '''Sets up self.info['Basefile'].'''
        ### Basefile and ResDir ###
        if self.info['Basefile'].lower() in ['','none']:
            self.info['Basefile'] = rje.baseFile(self.obj['SeqList'].info['Name'])
            if self.info['ResDir'].lower() not in ['','none']:
                rje.mkDir(self,self.info['ResDir'])
                self.info['Basefile'] = self.info['ResDir'] + self.dataset()
        ### Input ###
        if not self.info['Input']: self.info['Input'] = self.obj['SeqList'].info['Name']
#########################################################################################################################
    def seqBaseFile(self):  ### Returns the results directory basefile as specified by sequence alone (for SLiMSearch)
        '''Returns the results directory basefile as specified by sequence alone (e.g. for SLiMSearch)'''
        return self.info['ResDir'] + rje.baseFile(self.obj['SeqList'].info['Name'],strip_path=True)
#########################################################################################################################
    def setAlphabet(self):  ### Sets up self.list['Alphabet']
        '''Sets up self.list['Alphabet'].'''
        ### ~ [1] ~ Setup basic alphabet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.list['Alphabet']: self.list['Alphabet'] = string.split(string.join(self.list['Alphabet']).upper())
        else:
            if self.opt['DNA']: self.list['Alphabet'] = rje_seq.alph_dna[:-1]
            else: self.list['Alphabet'] = rje_seq.alph_protx[:-1]
        ### ~ [2] ~ Remove masked AAs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for aa in self.list['AAMask']:
            if aa in self.list['Alphabet']: self.list['Alphabet'].remove(aa)
#########################################################################################################################
    def setQRegion(self):   ### Sets up Query Region masking
        '''Sets up Query Region masking.'''
        if 'QRegion' not in self.list: self.list['QRegion'] = []
        try:
            qreg = []; pairs = True
            for x in self.list['QRegion']:
                r = string.atoi(x)
                qreg.append(r)
                pairs = not pairs
            if not pairs: qreg.append(-1)
            self.list['QRegion'] = qreg
        except:
            self.errorLog('Problem with QRegion masking %s - will not use' % self.list['QRegion'])
            self.list['QRegion'] = []
        if self.list['QRegion']: self.printLog('#QREG','%d Query region(s)' % int((len(self.list['QRegion'])+0.5)/2))
#########################################################################################################################
    def maskInput(self):    ### Masks input sequences, replacing masked regions with Xs
        '''Masks input sequences, replacing masked regions with Xs. Creates seq.info['PreMask' & 'MaskSeq']'''
        ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.maskPickle(load=True): return
        masked = False
        if not self.opt['Masking']:
            for seq in self.seqs(): seq.info['PreMask'] = seq.info['MaskSeq'] = seq.info['Sequence'][0:]
            self.opt['Masked'] = False
            return
        ## ~ [0a] ~ UniFake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.opt['FakeMask']:
            fake = unifake.UniFake(self.log,self.cmd_list)
            fake.setup()
            fake.obj['SeqList'] = self.obj['SeqList']
            fake.info['DatOut'] = '%s.dat' % self.info['Basefile']
            fake.uniFake(store=True)
            try: self.printLog('#FAKE','%d UniFake entries added for masking' % len(self.obj['SeqList'].obj['UniProt'].list['Entry']))
            except: self.errorLog('UniFake Masking error')
            self.wallTime()
        ### ~ [1] ~ Case Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for seq in self.seqs():
            seq.info['PreMask'] = seq.info['Sequence'][0:]
            if seq.maskCase(case=self.info['CaseMask'],mask='X',log=self.opt['LogMask']): masked = True
            self.wallTime()
        ### ~ [2] ~ Disorder (Inclusive Filtering) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['DisMask']:
            for seq in self.seqs():
                seq.maskDisorder(inverse=True,log=self.opt['LogMask'])
                masked = True
                self.wallTime()
        ### ~ [3] ~ Relative conservation masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['ConsMask']:
            slimcalc = rje_slimcalc.SLiMCalc(self.log,self.cmd_list+['conspec=','slimcalc=','usealn=T'])
            slimcalc.loadBLOSUM()
            cmx = 0
            for seq in self.seqs():
                (seq.info['MaskSeq'],seq.info['Sequence']) = (seq.info['Sequence'],seq.info['PreMask'])
                seq.deGap()
                relcon = slimcalc.relConListFromSeq(seq,window=30)
                newseq = []
                try:
                    for i in range(seq.seqLen()):
                        if relcon[i] < 0: newseq.append('X')
                        else: newseq.append(seq.info['MaskSeq'][i])
                    maskx = newseq.count('X') - seq.info['MaskSeq'].count('X'); cmx += maskx
                    seq.info['Sequence'] = string.join(newseq,'')
                    if maskx > 0 and self.opt['LogMask']: self.printLog('#MASK','Masked %s low relative conservation. (%d X added to %daa seq.)' % (seq.shortName(),maskx,seq.aaLen()))
                except:
                    self.errorLog('Problem with relative conservation masking for %s' % seq.shortName())
                    self.printLog('#REL','%daa (%d) => %d conscores' % (seq.seqLen(),seq.aaLen(),len(relcon)))
                    seq.info['Sequence'] = seq.info['MaskSeq']
                self.wallTime()
            self.printLog('#REL','%s aa masked using relative conservation' % rje.integerString(cmx))
        ### ~ [4] ~ UniProt Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.obj['SeqList'].obj['UniProt']:
            for entry in self.obj['SeqList'].obj['UniProt'].list['Entry']:
                if self.list['IMask']:
                    entry.maskFT(types=self.list['IMask'],inverse=True,mask='X',log=self.opt['LogMask'])
                    masked = True
                if self.list['FTMask']:
                    entry.maskFT(types=self.list['FTMask'],inverse=False,mask='X',log=self.opt['LogMask'])
                    masked = True
                self.wallTime()
        ### ~ [5] ~ Low Complexity Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if rje.matchExp('(\d+),(\d+)',self.info['CompMask']):
            (x,y) = rje.matchExp('(\d+),(\d+)',self.info['CompMask'])
            for seq in self.seqs(): seq.maskLowComplexity(lowfreq=int(x),winsize=int(y),log=self.opt['LogMask'])
            masked = True
            self.wallTime()
        elif self.info['CompMask'].lower()[:1] not in ['f','n']:
            self.log.errorLog('CompMask "%s" wrong format for complexity masking' % self.info['CompMask'],printerror=False)
        ### ~ [6] ~ N-terminal Met ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['MaskM']:
            mx = 0
            for seq in self.seqs():
                if seq.info['Sequence'][0] == 'M':
                    seq.info['Sequence'] = 'X' + seq.info['Sequence'][1:]
                    mx += 1
            self.printLog('#MASK','MetMask: %d N-terminal Ms masked' % mx)
            self.wallTime()
        ### ~ [7] ~ Position-specific masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.dict['MaskPos']:
            for pos in rje.sortKeys(self.dict['MaskPos']):
                try: (x,y) = (int(pos),rje.strList(self.dict['MaskPos'][pos]))
                except: self.log.errorLog('Problem with MaskPos entry "%s:%s"' % (pos,self.dict['MaskPos'].pop(pos)))
            for seq in self.seqs(): seq.maskPosAA(self.dict['MaskPos'],mask='X',log=self.opt['LogMask'])
            masked = True; self.wallTime()
        ### ~ [8] ~ Motif masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.info['MotifMask']:
            maskslims = rje_slimlist.SLiMList(self.log,self.cmd_list)
            maskslims.loadMotifs(self.info['MotifMask'])
            mx = 0
            for seq in self.seqs():
                sx = -seq.info['Sequence'].count('X')
                for slim in maskslims.slims():
                    for hit in slim.searchSequence(sequence=seq.info['PreMask']):    # {Pos,Variant,Match,ID,MisMatch}
                        #self.deBug(hit)
                        try:
                            r = hit['Pos'] - 1
                            for i in range(len(hit['Variant'])):
                                if hit['Variant'][i] != '.': seq.info['Sequence'] = rje.strSub(seq.info['Sequence'],r+i,r+i,'X')
                        except: self.errorLog('Problem with motifmask %s' % hit)
                sx += seq.info['Sequence'].count('X')
                if sx > 0:
                    mx += sx
                    if self.opt['LogMask']: self.printLog('#MASK','Motif-Masked %s. (%d X added to %daa seq.)' % (seq.shortName(),sx,seq.aaLen()))
                self.wallTime()
            self.printLog('#MASK','MotifMask: %d motif match residues masked' % mx)
        ### ~ [9] ~ AA Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.list['AAMask']:
            mx = 0
            for seq in self.seqs(): mx += seq.maskAA(self.list['AAMask'],log=self.opt['LogMask'])
            self.printLog('#MASK','AA-Masked %s: %s X added.' % (string.join(self.list['AAMask'],';'),rje.integerString(mx)))
            self.wallTime()
        ### ~ [10] ~ QRegion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:
            self.deBug(self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus'])
            self.deBug(self.list['QRegion'])
            self.deBug('Focus' in self.dict)
            self.deBug('Query' in self.dict['Focus'])
        except: pass
        if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']:
            mx = 0
            for seq in self.dict['Focus']['Query']: mx += seq.maskRegion(self.list['QRegion'],inverse=True)
            self.printLog('#MASK','QRegion %s: %s X added.' % (self.list['QRegion'],rje.integerString(mx)))
            self.wallTime()
        ### ~ [11] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for seq in self.seqs(): seq.info['MaskSeq'] = seq.info['Sequence'][0:]
        if masked:
            self.log.printLog('#MASK','%s: Masking of input sequences complete.' % self.dataset())
            self.opt['Masked'] = True
        #self.deBug(self.obj['SeqList'].seq[0].info)
        self.maskPickle()
#########################################################################################################################
    def maskPickle(self,load=False):    ### Loads or Saves post-masking pickle
        '''Loads or Saves post-masking pickle.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['MaskPickle']: return False
            maskpickle = '%s.mask.%s' % (rje.baseFile(self.info['Input'],strip_path=True),self.maskText(freq=False))  ## Pickle name ##
            
            ### ~ [1] Load existing pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if load:
                newme = None    ## New object, loaded from Pickle
                ## Check for file and load ##
                for pdir in [self.info['ResDir'],self.info['BuildPath']]:
                    pfile = '%s%s' % (pdir,maskpickle)
                    newme = self.unpickleMe(basefile=pfile,process=False)
                    if newme: break
                if not newme: return None
                ## Check for masking differences ##
                changes = []
                for var in ['CompMask','CaseMask','MotifMask']:         # Info
                    if self.info[var] != newme.info[var]: changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['Masking','DisMask','MaskM','ConsMask']:   # Opt
                    if self.opt[var] != newme.opt[var]: changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['FTMask','IMask']:                          # List   
                    slist = self.list[var][0:]
                    nlist = newme.list[var][0:]
                    slist.sort()
                    nlist.sort()
                    if slist != nlist:
                        changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                ## Recreate or use pickle but add new commands ##
                if changes and (self.stat['Interactive'] < 0 or rje.yesNo('%d SLiMBuild parameter mismatches with pickle. Create new pickle?' % len(changes))):
                    self.printLog('#PICKLE','Parameters changed. Making new pickle.')
                    return None
                newme.cmd_list = self.cmd_list
                newme.setInfo(self.info)
                newme.setStat(self.stat)
                self.opt['Masked'] = newme.opt['Masked']
                newme.setOpt(self.opt)
                newme.info['ResFile'] = self.info['ResFile']
                newme.info['ResDir'] = self.info['ResDir']
                newme.info['BuildPath'] = self.info['BuildPath']
                newme.stat['StartTime'] = self.stat['StartTime']
                newme.setLog(self.log)
                for o in self.obj:
                    if o != 'SeqList': newme.obj[o] = self.obj[o]
                self.list['UP'] = newme.list['UP']
                self.dict['AAFreq'] = newme.dict['AAFreq']
                self.dict['MST'] = newme.dict['MST']
                newme.list = self.list
                newme.dict = self.dict
                self.replaceMe(newme)
                #self.deBug(self.obj['SeqList'].seq[0].info)
                return True
                
            ### ~ [2] Save new pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.maskPickleMe(basefile='%s%s' % (self.info['ResDir'],maskpickle),gzip=True,replace=True)
            return None
        except: self.errorLog('Problem with %s.maskPickle()' % self.prog()); return False
#########################################################################################################################
    def maskPickleMe(self,basefile=None,gzip=True,replace=True):   ### Saves self object to pickle and zips
        '''
        Saves self object to pickle and zips.
        >> basefile:str [None] = if none, will use self.info['Basefile']
        >> gzip:bool [True] = whether to GZIP (win32=F only)
        >> replace:bool [True] = whether to replace existing Pickle
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Pickle']: self.printLog('#PICKLE','No pickling. (pickle=F)'); return False
            if not basefile or basefile.lower() == 'none': basefile = self.info['Basefile']
            pfile = '%s.pickle' % basefile
            if self.opt['Webserver']: pfile = '%s%s' % (self.info['BuildPath'],os.path.basename(pfile))
            if not replace and (os.path.exists(pfile) or os.path.exists('%s.gz' % pfile)):
                self.printLog('#PICKLE','No pickling - pickle already exists (no replace)!'); return False
            ### ~ [2] ~ Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SAVE','Attempting to save to %s.' % pfile,log=False)
            pickle.dump(self,open(pfile,'w'))
            self.printLog('#SAVE','Intermediate saved as %s (Python pickle).' % pfile)
            ### ~ [3] ~ GZip and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Win32'] and gzip:
                try:
                    if os.path.exists('%s.gz' % pfile): os.unlink('%s.gz' % pfile)
                    os.system('gzip %s' % pfile)
                    self.printLog('#GZIP','%s zipped.' % pfile)
                except: self.errorLog('Cannot gzip %s' % pfile)
            return True
        except: self.errorLog('Problem during %s.pickleMe()' % self); return False
#########################################################################################################################
    def makeUPC(self):  ### Generates UP Clusters from self.obj['SeqList'] using BLAST
        '''Generates UP Clusters from self.obj['SeqList'] using BLAST.'''
        try:### ~ [1] Check whether BLAST and MST is necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] No evolutionary filtering (not a good idea generally!) ~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.opt['EFilter']:
                self.list['UP'] = []
                for seq in self.obj['SeqList'].seq: self.list['UP'].append((seq,))
                self.makeMST()
                self.printLog('#EVOL','WARNING! No evolutionary filtering (efilter=F)')
                return True
            ## ~ [1b] Direct reading of *.upc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.opt['Force'] and self.readUPC(): return True
            self.list['UP'] = []
            self.dict['MST'] = {}            
            ### ~ [2] Setup BLAST etc - necessary for UPC and/or Focal Group MST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gentxt = ''
            upcbase = '%s%s' % (self.info['ResDir'],rje.baseFile(self.info['Input'],strip_path=True))
            self.obj['SeqList'].info['Name'] = seqfile = '%s.slimdb' % upcbase
            self.obj['SeqList'].saveFasta()
            seqdict = self.obj['SeqList'].seqNameDic()
            gablam = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            ## ~ [2a] Read distance matrix if present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            readdis = False
            gablamdis = False
            mfiles = ['%s%s.dis.tdt' % (self.info['ResDir'],rje.baseFile(self.info['Input'],strip_path=True)),
                      '%s%s.dis.tdt' % (self.info['BuildPath'],rje.baseFile(self.info['Input'],strip_path=True))]
            if self.opt['Force']: mfiles = []
            if self.info['AltDis'].lower() not in ['','none']: mfiles.append(self.info['AltDis'])
            if self.info['GablamDis'].lower() not in ['','none']: mfiles.append(self.info['GablamDis'])
            for dismatrix in mfiles:
                if os.path.exists(dismatrix):
                    if dismatrix == self.info['GablamDis']:
                        gablam.loadFromDataTable(dismatrix)
                        gentxt = 'GABLAM DisMat %s' % dismatrix
                    else:
                        gablam.loadMatrix(dismatrix,checksym=False,default=1.0)
                        gentxt = 'DisMat %s' % dismatrix
                    self.log.printLog('#UPC','Loaded distance matrix from %s' % dismatrix)
                    if dismatrix == self.info['GablamDis']: gablamdis = True
                    break
            if gablam.dict['Matrix']:   # Something read... convert
                updict = {}
                readdis = True
                newmatrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
                (cx,ctot) = (0.0,self.seqNum())
                for seq1 in self.seqs():
                    self.log.printLog('\r#UPC','Converting distance matrix: %.1f%%' % (cx/ctot),newline=False,log=False)
                    cx += 100.0
                    updict[seq1] = []
                    if seq1.shortName() not in gablam.dict['Matrix']: continue
                    for upseq in gablam.dict['Matrix'][seq1.shortName()]:   #x#for seq2 in self.seqs():
                        if upseq not in seqdict: continue
                        seq2 = seqdict[upseq]
                        try: newmatrix.addDis(seq1,seq2,gablam.dict['Matrix'][seq1.shortName()][seq2.shortName()])
                        except:
                            if not gablamdis: 
                                self.errorLog('Error converting distance matrix')
                                readdis = False
                                break
                        if newmatrix.getDis(seq1,seq2,default=1.0) < 1 or newmatrix.getDis(seq2,seq1,default=1.0) < 1: updict[seq1].append(seq2)
                    if not readdis:
                        del newmatrix
                        break
                self.log.printLog('\r#UPC','Convertion of distance matrix complete.',log=False)
            if readdis: gablam.dict['Matrix'] = newmatrix.dict['Matrix']
            else:
            ## ~ [2b] BLAST generate matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                blast = rje_blast.BLASTRun(self.log,['blastcf=F','blastf=T','blaste=1e-4']+self.cmd_list)
                blast.setInfo({'InFile':seqfile,'DBase':seqfile,'Name':'%s.self.blast' % upcbase})
                blast.setStat({'OneLine':self.seqNum(),'HitAln':self.seqNum()})
                if self.opt['DNA']: blast.info['Type'] = 'blastn'
                blast.formatDB(fasfile=seqfile,force=self.opt['Force'],protein=not self.opt['DNA'])
                if not rje_blast.checkForDB(dbfile=seqfile,checkage=False,log=self.log,protein=not self.opt['DNA']):
                    self.errorLog('FormatDB failed for unknown reasons.',printerror=False)
                    raise IOError
                ### Perform and read BLAST ###
                blast.blast(cleandb=True,use_existing=not self.opt['Force'])
                if not blast.readBLAST(gablam=True,unlink=not self.extras(2)):
                    self.log.errorLog('Major problem with BLAST for unknown reasons.')
                    raise IOError
                gentxt = 'blaste=%s, blastcf=%s, blastf=%s' % (rje_slim.expectString(blast.stat['E-Value']),str(blast.opt['Composition Statistics'])[0],str(blast.opt['Complexity Filter'])[0])
                ### Make GABLAM Distance Matrix from BLAST ###
                updict = {}
                for seq in self.seqs(): gablam.addDis(seq,seq,0.0)      # Always self-distance of zero, even without BLAST hits
                for search in blast.search:
                    seq = seqdict[search.info['Name']]
                    updict[seq] = []
                    for hit in search.hit:
                        if seqdict.has_key(hit.info['Name']):
                            updict[seq].append(seqdict[hit.info['Name']])
                            if seqdict[hit.info['Name']] != seq:
                                gablam.addDis(seq,seqdict[hit.info['Name']],1 - (hit.dict['GABLAM']['Query']['GABLAM ID'] / float(seq.aaLen())))
                        else: self.errorLog('No sequence for hit "%s"!' % hit.info['Name'],printerror=False)
            ## ~ [2c] Check and enforce symmetry of UP information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gablam.forceSymmetry(method='min',missing=1.0)
            for seq in rje.sortKeys(updict):
                for hom in updict[seq][0:]:
                    while updict[seq].count(hom) > 1: updict[seq].remove(hom) 
                    if seq not in updict[hom]: updict[hom].append(seq)

            ##!# New DomTable supplement of GABLAM #!##
            if 'DomTable' not in self.info: self.info['DomTable'] = ''
            if rje.checkForFile(self.info['DomTable']):
                domdata = {}
                seqdom = rje.dataDict(self,self.info['DomTable'],['Name'],['Type'],lists=True)
                self.printLog('#DOM','Using domain table for additional UPC grouping.')
                for seq in self.seqs():
                    if seq.shortName() in seqdom:
                        for dom in seqdom[seq.shortName()]['Type']:
                            if dom not in domdata: domdata[dom] = []
                            domdata[dom].append(seq)
                #self.deBug(domdata)
                for dom in domdata:
                    for seq1 in domdata[dom]:
                        for seq2 in domdata[dom]:
                            if seq1 != seq2 and seq1 not in updict[seq2]: updict[seq2].append(seq1)
            elif self.info['DomTable'].lower() not in ['','none']: self.errorLog('Domain table "%s" not found' % self.info['DomTable'],printerror=False)

            ##!# New HomCut reduction of sequences based on BLAST #!##
            if self.stat['HomCut']:
                (px,prex) = (0.0,self.seqNum())
                ## Remove hub proteins ##
                for seqname in rje.sortKeys(seqdict):
                    self.log.printLog('\r#HOM','Checking HomCut (%d) %.1f%%' % (self.stat['HomCut'],px/prex),log=False,newline=False)
                    px += 100.0
                    seq = seqdict[seqname]
                    if len(updict[seq]) > self.stat['HomCut']:
                        updict.pop(seq)
                        seqdict.pop(seqname)
                        self.obj['SeqList'].removeSeq(text='Exceeds HomCut (%d)' % self.stat['HomCut'],seq=seq)
                self.log.printLog('\r#HOM','Checked HomCut (%d). %s of %s seqs remain.' % (self.stat['HomCut'],rje.integerString(self.seqNum()),rje.integerString(prex)))
                gentxt += '; HomCut=%d' % self.stat['HomCut']
                ## Clean up ##
                if self.seqNum() < prex:
                    (px,ptot) = (0.0,self.seqNum())
                    for seq in updict:
                        self.log.printLog('\r#HOM','Cleaning up after HomCut %.1f%%' % (px/ptot),log=False,newline=False)
                        px += 100.0
                        for hit in updict[seq][0:]:
                            if hit not in self.seqs(): updict[seq].remove(hit)
                    self.verbose(0,5,'\r',0)
                    self.obj['SeqList'].saveFasta()
                        
            ### Make UP Clusters ##
            self.printLog('#UP','UP clusters',log=False,newline=False)
            seqdict = self.obj['SeqList'].seqNameDic()
            sortedseq = rje.sortKeys(seqdict)
            while updict:
                sname = sortedseq.pop(0); seq = seqdict[sname]
                while seq not in updict: sname = sortedseq.pop(0); seq = seqdict[sname]
                self.list['UP'].append((seq,))
                uplen = 0
                while uplen < len(self.list['UP'][-1]):
                    self.printLog('\r#UP','%d UP clusters (%d seq remaining)     ' % (self.UPNum(),len(updict)),log=False,newline=False)
                    uplen = len(self.list['UP'][-1])
                    for seq in self.list['UP'][-1]:
                        if updict.has_key(seq):
                            for p in updict.pop(seq):
                                if p not in self.list['UP'][-1]:
                                    self.list['UP'][-1] += (p,)
            self.printLog('\r#UP','%d UP clusters generated from %s sequences.' % (self.UPNum(),rje.integerString(self.seqNum())),log=False)

            ### Make MST ###
            if gablam and gablamdis: self.meanDisMST(gablam)
            else: self.makeMST(gablam)

            ### Output *.upc ###
            UFILE = open('%s.upc' % upcbase,'w')
            UFILE.write('#%s# %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.stat['MST'],gentxt))
            rje.writeDelimit(UFILE,['UP','N','MST','Seqs'],'\t')
            for u in range(self.UPNum()):
                upc = self.list['UP'][u]
                seqs = upc[0].shortName()
                for seq in upc[1:]: seqs += ' %s' % seq.shortName()
                rje.writeDelimit(UFILE,[u+1,len(upc),'%.3f' % (len(upc)*self.dict['MST'][upc]),seqs],'\t')
            UFILE.close()
            self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.stat['MST'],gentxt))

            ### Output DisMatrix ###
            gablam.info['Name'] = '%s GABLAM' % self.dataset()
            if self.stat['MaxUPC'] > 0 and self.stat['MaxUPC'] < self.UPNum():
                self.log.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). No *.dis.tdt.' % (self.UPNum(),self.stat['MaxUPC']))
            else:
                gablam.saveMatrix(self.seqs(),filename='%s.dis.tdt' % upcbase,delimit='\t',format='text',log=True)
                if self.extras(2):
                    gablam.saveMatrix(self.seqs(),filename='%s.phydis.txt' % upcbase,format='phylip',log=True)
            return True                
                
        except: self.log.errorLog('Fatal Error in %s.makeUPC().' % self.prog())
        return False
#########################################################################################################################
    def readUPC(self):  ### Looks for UPC file and loads details from file.
        '''Looks for UPC file and loads details from file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqbase = self.info['ResDir'] + rje.baseFile(self.obj['SeqList'].info['Name'],strip_path=True)
            self.obj['SeqList'].info['Name'] = seqfile = '%s.slimdb' % self.info['Basefile']    # Replace with seqbase
            self.list['UP'] = []
            self.dict['MST'] = {}
            seqdict = self.obj['SeqList'].seqNameDic()
            ulines = []
            bfiles = [seqbase,'%s%s' % (self.info['BuildPath'],os.path.split(seqbase)[1]),self.info['Basefile'],'%s%s' % (self.info['BuildPath'],self.dataset())]
            if bfiles[0] == bfiles[1]: bfiles.pop(1)
            for pfile in bfiles:
                ufile = '%s.upc' % pfile
                if os.path.exists(ufile): return self.loadUPCFromFile(ufile)
                self.printLog('#UPC','UPC file "%s" not found.' % ufile)
            return False
        except ValueError: self.errorLog('Problem with *.upc formatting in %s.readUPC().' % self.prog())
        except: self.errorLog('Error in %s.readUPC().' % self.prog())
        return False
#########################################################################################################################
    def loadUPCFromFile(self,ufile):    ### Loads UPC details from specific file
        '''Loads UPC details from specific file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(ufile):
                self.printLog('#UPC','Loading UPC from "%s".' % ufile)
                ulines = self.loadFromFile(ufile,chomplines=True)
            else:
                self.printLog('#UPC','UPC file "%s" not found.' % ufile)
                return False
            seqdict = self.obj['SeqList'].seqNameDic()
            ### ~ [1] ~ Read UPC File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['UP'] = []
            self.dict['MST'] = {}
            ## ~ [1a] ~ Dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('\r#UPC',ulines[0])
            data = rje.matchExp('#(\S.*)# (\d+) Seq; (\d+) UPC; (\S+) MST',ulines[0])
            base = data[0]
            seqx = int(data[1])
            if seqx != self.seqNum(): self.printLog('#ERR','Wrong number of sequences in UPC file!'); return False
            upx = int(data[2])
            self.stat['MST'] = float(data[3])
            ## ~ [1b] ~ UPCs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for u in range(upx):
                data = string.split(ulines[u+2],'\t')
                if int(data[0]) != (u+1): raise ValueError
                useqs = []
                for seq in string.split(data[-1]): useqs.append(seqdict[seq])
                if len(useqs) != int(data[1]): raise ValueError
                upc = (useqs[0],)
                for seq in useqs[1:]: upc += (seq,)
                if len(upc) != int(data[1]): raise ValueError
                seqx -= len(upc)
                self.dict['MST'][upc] = float(data[2]) / len(upc)
                self.list['UP'].append(upc)
            if upx != self.UPNum(): self.printLog('#ERR','Wrong number of UPC in UPC file!'); return False
            if seqx != 0: self.printLog('#ERR','Wrong number of sequences read from UPC file!'); return False
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.stat['MST']))
            return True        
        except ValueError: self.errorLog('Problem with *.upc formatting or content in %s.loadUPCFromFile().' % self.prog())
        except: self.errorLog('Error in %s.loadUPCFromFile().' % self.prog())
        return False
#########################################################################################################################
    def makeMST(self,gablam=None):   ### Makes UPC dictionary from GABLAM
        '''Makes UPC dictionary from GABLAM.'''
        self.dict['MST'] = {}
        if gablam:
            gablam.checkSymmetry(force=True)
            self.stat['MST'] = gablam.MST()
            for upc in self.list['UP']:
                uplist = []
                for seq in upc: uplist.append(seq)
                self.dict['MST'][upc] = gablam.MST(uplist) / len(upc)   # Mean Total * MST Value
        else:
            self.stat['MST'] = self.seqNum()
            for upc in self.list['UP']: self.dict['MST'][upc] = 1.0
#########################################################################################################################
    def meanDisMST(self,gablam=None):   ### Makes UPC dictionary from GABLAM using mean distance instead of MST
        '''Makes UPC dictionary from GABLAM using mean distance instead of MST.'''
        self.dict['MST'] = {}
        self.stat['MST'] = 0.0
        self.log.printLog('#SYM','Checking symmetry of GABLAM DisMatrix')
        gablam.checkSymmetry(force=True)
        (ux,upx) = (0.0,100.0/len(self.list['UP']))
        for upc in self.list['UP']:
            self.log.printLog('\r#MST','Calculating MeanDis MST replacement %.1f%%' % ux,newline=False,log=False)
            ux += upx
            uplist = []
            if len(upc) == 1:
                self.dict['MST'][upc] = 1.0
                self.stat['MST'] += 1.0
                continue
            self.dict['MST'][upc] = float((len(upc) * (len(upc)-1)))
            for s1 in upc:
                if s1 not in gablam.dict['Matrix']: continue
                for s2 in gablam.dict['Matrix'][s1]:
                    if s1 == s2 or s2 not in upc: continue
                    self.dict['MST'][upc] -= gablam.dict['Matrix'][s1][s2]
            self.dict['MST'][upc] = self.dict['MST'][upc] / (len(upc) * (len(upc)-1))   # Mean Total * MST Value
            self.stat['MST'] += self.dict['MST'][upc] * len(upc)
        self.log.printLog('\r#MST','Calculation of MeanDis MST replacement complete.')
#########################################################################################################################
    ### <5> ### SLiM Calculation method                                                                                 #
#########################################################################################################################
    def addSLiMToList(self,slim):   ### Add slims from self.dict['Slim'] to self.obj['SlimList']
        '''Add slims from self.dict['Slim'] to self.obj['SlimList'].'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            Motif = rje.dictValues(self.dict['Slim'][slim],'SLiM','SLiM')
            if Motif: return Motif
            pattern = patternFromCode(slim)
            Motif = self.obj['SlimList']._addMotif(slim,pattern,reverse=False,check=False,logrem=True)
            self.dict['Slim'][slim]['SLiM'] = Motif

            ### ~ [2] Add Calculate Statistics and make MotifOcc objects ###
            for (Seq,pos) in self.dict['Slim'][slim]['Occ'][0:]:
                if Seq not in Motif.dict['Occ']: Motif.dict['Occ'][Seq] = []
                if pos < 0: pos = 0     ## Re-find motif in sequence to get proper object data ##
                sequence = Seq.info['PreMask'][pos:pos+Motif.slimLen()]     #x#stat['FullLength']]
                #x#self.deBug('%s %s : %s\n%s' % (pos, slim, sequence,Seq.info['PreMask'][:pos+Motif.slimLen()]))
                for occ in Motif.searchSequence(sequence=sequence):    # {Pos,Variant,Match,ID,MisMatch}
                    try:
                        occ['Start_Pos'] = occ['Pos'] + pos
                        if slim[0] == '^': occ['Start_Pos'] += 1
                        occ['End_Pos'] = occ['Start_Pos'] + len(occ['Match']) - 1
                        occ['Expect'] = self.dict['Slim'][slim]['ExpUP']
                    except:
                        print Seq.info['Sequence']
                        print Seq.info['PreMask']
                        print slim, pos, Seq.shortName(), sequence[:10], Motif.dict['Search']
                        occ = {'Match':'!ERR!','SearchVar':'!ERR!','Variant':'!ERR!'}
                    occ['Pos'] = pos+1
                    occ['Prot_Len'] = Seq.aaLen()
                    occ['Seq'] = Seq
                    occ['Motif'] = Motif
                    #self.deBug(occ)
                    Motif.dict['Occ'][Seq].append(occ)
            return Motif
        except:
            self.errorLog('Problem with %s.addSLiMToList(%s)' % (self.prog(),slim))
            return None
#########################################################################################################################
    def OLDaddSLiMToList(self,slim):   ### Add slims from self.dict['Slim'] to self.obj['SlimList']
        '''Add slims from self.dict['Slim'] to self.obj['SlimList'].'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            Motif = rje.dictValues(self.dict['Slim'][slim],'SLiM','SLiM')
            if Motif: return Motif
            pattern = patternFromCode(slim)
            Motif = self.obj['SlimList']._addMotif(slim,pattern,reverse=False,check=False,logrem=True)
            self.dict['Slim'][slim]['SLiM'] = Motif

            ### ~ [2] Add Calculate Statistics and make MotifOcc objects ###
            for (Seq,pos) in self.dict['Slim'][slim]['Occ'][0:]:
                if Seq not in Motif.dict['Occ']: Motif.dict['Occ'][Seq] = []
                if pos < 0: pos = 0     ## Re-find motif in sequence to get proper object data ##
                sequence = Seq.info['PreMask'][pos:pos+Motif.slimLen()]     #x#stat['FullLength']]
                #x#self.deBug('%s %s : %s\n%s' % (pos, slim, sequence,Seq.info['PreMask'][:pos+Motif.slimLen()]))
                try:
                    occ = Motif.searchSequence(sequence=sequence)[0]    # {Pos,Variant,Match,ID,MisMatch}
                    occ['Start_Pos'] = occ['Pos'] + pos
                    if slim[0] == '^': occ['Start_Pos'] += 1
                    occ['End_Pos'] = occ['Start_Pos'] + len(occ['Match']) - 1
                    occ['Expect'] = self.dict['Slim'][slim]['ExpUP']
                except:
                    print Seq.info['Sequence']
                    print Seq.info['PreMask']
                    print slim, pos, Seq.shortName(), sequence[:10], Motif.dict['Search']
                    occ = {'Match':'!ERR!','SearchVar':'!ERR!','Variant':'!ERR!'}
                occ['Pos'] = pos+1
                occ['Prot_Len'] = Seq.aaLen()
                occ['Seq'] = Seq
                occ['Motif'] = Motif
                self.deBug(occ)
                Motif.dict['Occ'][Seq].append(occ)
            return Motif
        except:
            self.errorLog('Problem with %s.addSLiMToList(%s)' % (self.prog(),slim))
            return None
#########################################################################################################################
    def calculateSLiMOccStats(self):   ### Makes entries to SLiMList object and calculates attributes with slimcalc
        '''Makes entries to SLiMList object and calculates attributes with slimcalc.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for slim in self.dict['Slim']: self.addSLiMToList(slim)
            self.obj['SlimList'].calculateOccAttributes(silent=False)
            self.opt['OccStatsCalculated'] = True
        except: self.log.errorLog('Problem with %s.calculateSLiMOccStats()' % self.prog())
#########################################################################################################################
    ### <6> ### SLiMChance Probability Methods                                                                          #
#########################################################################################################################
    def makeAAFreq(self):   ### Makes an initial AAFreq dictionary containing AA counts only (including Xs)
        '''Makes an initial AAFreq dictionary containing AA counts only (including Xs).'''
        try:### ~ [1] ~ Load AA Frequencies for file and use for whole dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['AAFreq'].lower() not in ['','none']:
                self.opt['MaskFreq'] = False     # Over-riding AA Frequencies invalidates premasking
                self.opt['SmearFreq'] = False   # Amino acids are already constant across UPC
                ## Read Freqs ##
                aafile = self.info['AAFreq']
                if os.path.exists(aafile) and open(aafile,'r').read()[:1] == '>':   # Fasta
                    aafreq = self.obj['SeqList'].aaFreq(fromfile=aafile,total=True)    # Returns proportions and total
                else: aafreq = self.obj['SeqList'].aaFreq(loadfile=aafile,total=True)    # Returns proportions and total
                ## Adjust to Counts ##
                for aa in aafreq.keys():
                    if aa != 'Total': aafreq[aa] = int(aafreq[aa]*aafreq['Total']+0.5)
                ## Copy to seqs and UPCs ##
                for dkey in ['Dataset'] + self.list['UP']: self.dict['AAFreq'][dkey] = copy.copy(aafreq)
                self.dict['AAFreq']['Dataset']['Total'] = self.obj['SeqList'].aaTotal(nonx=True)  ### Returns total number of AA in dataset
                self.printLog('#AAFREQ','Using aa frequencies from %s (no mask/smear freq)' % self.info['AAFreq'])
            ### ~ [2] ~ Calculate from Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                self.dict['AAFreq']['Dataset'] = {}
                for upc in self.list['UP']:
                    self.dict['AAFreq'][upc] = {}
                    for seq in upc:
                        seq.aaFreq(aafreq=self.dict['AAFreq'][upc])
                        seq.aaFreq(aafreq=self.dict['AAFreq']['Dataset'])
                ## Add Totals to AAFreq dict ## 
                for dkey in ['Dataset'] + self.list['UP']:    #!# Don't need totals? #!#
                    self.dict['AAFreq'][dkey]['Total'] = sum(self.dict['AAFreq'][dkey].values())
            ### ~ [3] ~ Convert Dataset into Freqs - not bothered with masking for this ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.dictFreq(self.dict['AAFreq']['Dataset'])
        except: self.errorLog('Error in SLiMFinder.makeAAFreq().');  raise
#########################################################################################################################
    def smearAAFreq(self,update=True):  ### Equalises AAFreq across UPC. Leaves Totals unchanged.
        '''Equalises AAFreq across UPC. Leaves Totals unchanged. Updates or returns smearfreq.'''
        try:
            ###~Setup dictionary containing mean AAFreq across UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            smearfreq = {}
            letters = {False:'ACDEFGHIKLMNPQRSTVWY^$',True:'ACGT^$'}[self.opt['DNA']]
            for aa in letters:
                smearfreq[aa] = 0.0
                for upc in self.list['UP']: smearfreq[aa] += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)
            rje.dictFreq(smearfreq,total=False)
            if not update: return smearfreq
            ###~Copy frequencies to all AAFreq dictionaries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            for grp in self.dict['AAFreq']:
                for aa in letters: self.dict['AAFreq'][grp][aa] = smearfreq[aa]
            if self.opt['DNA']: self.printLog('#NTFREQ','NT frequencies "smeared" over %d UPC' % self.UPNum())
            else: self.printLog('#AAFREQ','AA frequencies "smeared" over %d UPC' % self.UPNum())
        except:
            self.errorLog('Major problem during %s.smearAAFreq()' % self.prog())
            raise
#########################################################################################################################
    ### <7> ### Results Output Methods                                                                                  #
#########################################################################################################################
    def extras(self,level=1):     ### Returns whether extras settings permit output at this level
        '''Returns whether extras settings permit output at this level.'''
        if 'Extras' in self.stat: return self.stat['Extras'] >= level
        return self.opt['Extras']
#########################################################################################################################
    def maskText(self,joiner='',freq=True): ### Returns masking text 
        '''Returns masking text.'''
        if self.info['MaskText'].lower() not in ['','none']: return self.info['MaskText']
        if not 'MaskText' in self.list: 
            masking = []
            if self.opt['Masking']:
                if self.opt['ConsMask']: masking.append('Cons')
                if self.opt['DisMask']: masking.append('Dis')
                if rje.matchExp('(\d+),(\d+)',self.info['CompMask']) and self.info['CompMask'] != '0,0':
                    masking.append('Comp-%s-%s' % rje.matchExp('(\d+),(\d+)',self.info['CompMask']))
                if self.list['FTMask']:
                    if self.opt['Webserver']:
                        ftxt = []
                        for f in self.list['FTMask']: ftxt.append(f[:1])
                        ftxt.sort()
                        masking.append('FT-%s' % string.join(ftxt,''))
                    else: masking.append('FT')
                if self.list['IMask']:
                    if self.opt['Webserver']:
                        ftxt = []
                        for f in self.list['IMask']: ftxt.append(f[:1])
                        ftxt.sort()
                        masking.append('Inc-%s' % string.join(ftxt,''))
                    else: masking.append('Inc')
                if self.info['MotifMask']: masking.append('Mot')
                if self.list['AAMask']: masking.append('AA')
                if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']: masking.append('QReg')
            if not masking: masking = ['NoMask']
            self.list['MaskText'] = masking
        if freq and self.opt['MaskFreq'] and 'Freq' not in self.list['MaskText']: self.list['MaskText'].insert(0,'Freq')
        if not freq and 'Freq' in self.list['MaskText']: self.list['MaskText'].remove('Freq')
        return string.join(self.list['MaskText'],joiner)
#########################################################################################################################
    def tidyMotifObjects(self): ### This should not be necessary but somehow is!
        '''This should not be necessary but somehow is!'''
        return #!# Try without this! #!#
        ### Tidy MotifList Object ###
        for Motif in self.obj['MotifList'].motifs()[0:]:
            slim = Motif.slimCode() 
            if not slim or slim not in self.list['SigSlim']: self.obj['MotifList'].removeMotif(Motif)
        ### Check completeness of MotifObject ###
        for slim in self.list['SigSlim']:   #!# Add SigNum? #!#
            pattern = patternFromCode(slim)
            Motif = self.obj['MotifList'].mapPattern(pattern,update=False)
            if not Motif in self.obj['MotifList'].motifs():
                self.log.printLog('#MOT','Motif "%s" missing from Motif Occurrence dictionary.' % pattern)
                try:
                    Motif = self.obj['MotifList'].mapPattern(pattern,update=True)
                    Motif = self.motifOccStats(slim)
                    Motif.stat['OccNum'] = self.slimOccNum(slim)
                    Motif.stat['OccSeq'] = self.slimUP(slim)
                    Motif.dict['Expect'] = {self.info['Basefile']:self.dict['Slim'][slim]['ExpUP']}
                    self.obj['MotifList'].combMotifOccStats(statlist=self.list['OccStats'],revlist=['Hyd'],log=False,motiflist=[Motif])
                    self.log.printLog('#MOT','Successfully added Motif Object for "%s" (%s)!' % (slim,pattern))
                except: self.errorLog('Cannot add/process Motif Object for "%s" (%s)!' % (slim,pattern))
#########################################################################################################################
    def abNprob(self,a,b,N,overlap=0,dirn='less'):     ### Returns probabilities of overlap in a given b/N and in b given a/N
        '''Returns probabilities of no overlap in a given b/N and in b given a/N.'''
        p = 0.0
        if dirn == 'less':
            for k in range(overlap+1):
                p += rje.binomial(k,a,float(b)/N,exact=True,callobj=self)
                p += rje.binomial(k,b,float(a)/N,exact=True,callobj=self)
        else: p += rje.binomial(overlap,b,float(a)/N,callobj=self) + rje.binomial(overlap,a,float(b)/N,callobj=self)
        return p / 2.0
#########################################################################################################################
    def myPickle(self):  ### Returns pickle identifier, also used for Outputs "Build" column (self.info['Build'])
        '''Returns pickle identifier, also used for Outputs "Build" column.'''
        ## Pickle name ##
        return '%s.%s' % (self.info['Build'],self.maskText())
#########################################################################################################################
    def tarZipSaveSpace(self):  ### Tars and Zips output and/or deletes extra files as appropriate
        '''Tars and Zips output and/or deletes extra files as appropriate.'''
        try:### ~ [0] Setup file lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            targz = '%s.%s.tar.gz' % (self.info['Basefile'],self.myPickle())
            ## ~ [0a] Pickle for this run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mypickle = '%s.%s.pickle' % (self.info['Basefile'],self.myPickle())
            gzpickle = '%s.gz' % mypickle
            if rje.isYounger(mypickle,gzpickle) == mypickle: os.unlink(gzpickle)
            if rje.isYounger(gzpickle,mypickle) == gzpickle: os.unlink(mypickle)
            if os.path.exists(mypickle) and not self.opt['Win32']:
                try:
                    os.system('gzip %s' % (mypickle))
                    self.printLog('#GZIP','%s %s zipped.' % (self.prog(),mypickle))
                except: self.errorLog('Cannot gzip %s' % (mypickle))
            ## ~ [0b] All dataset files in directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            files = glob.glob('%s*' % self.info['Basefile'])   
            ### ~ [1] Tar and Zip Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['TarGZ']:
                ## ~ [1a] Check Windows ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.opt['Win32']:
                    self.printLog('#TAR','Sorry! TarGZ option not available for Windows yet!')
                    if self.stat['SaveSpace'] > 0: self.printLog('#DEL','No files deleted to save space. Control output with extras=X.')
                    return
                ## ~ [1b] Delete existing archive ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if os.path.exists(targz):
                    os.unlink(targz)
                    self.printLog('#TAR','Deleted old TAR archive.')
                ## ~ [1c] Make archive ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #tcmd = 'tar -czf %s %s*' % (targz,self.info['Basefile'])
                tarfiles = []
                for tfile in files[0:]:
                    if tfile[-6:] == 'tar.gz': continue
                    if tfile.find('.pickle') > 0: continue
                    tarfiles.append(tfile)
                #for tfile in glob.glob('%s*tar.gz' % self.info['Basefile']):
                #    if tfile in tarfiles: tarfiles.remove(tfile)
                #for pfile in glob.glob('%s*pickle*' % self.info['Basefile']):
                #    if pfile in tarfiles: tarfiles.remove(pfile)
                if self.stat['SaveSpace'] > 2 and os.path.exists(mypickle): tarfiles.append(mypickle)
                if self.stat['SaveSpace'] > 2 and os.path.exists(gzpickle): tarfiles.append(gzpickle)
                if tarfiles:
                    tcmd = 'tar -czf %s %s' % (targz,string.join(tarfiles))
                    self.printLog('#TAR',tcmd)
                    if os.system(tcmd):
                        self.printLog('#ERR','Problem executing TarGZ option for %s!' % self.dataset())
                        if self.stat['SaveSpace'] > 0: self.printLog('#DEL','No files deleted to save space.')
                        return
            ### ~ [2] Delete files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['SaveSpace'] > 0:
                keeplist = ['tar.gz']
                if self.stat['SaveSpace'] < 3:  keeplist += ['upc','pickle*']
                if self.stat['SaveSpace'] < 2:  keeplist += ['occ.csv']
                ## ~ [2a] Setup File List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for keepext in keeplist:
                    keepers = glob.glob('%s*%s' % (self.info['Basefile'],keepext))
                    for k in keepers:
                        if k in files: files.remove(k)
                self.printLog('#DEL','Deleting of %d files/dirs to save space...' % len(files),newline=False,log=False)
                ## ~ [2b] Delete Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for file in files:
                    if os.path.isdir(file): rje.deleteDir(self,file,contentsonly=False,confirm=False,report=False)
                    else: os.remove(file)
                self.printLog('\r#DEL','Deletion of %d files/dirs to save space complete.' % len(files))
        except: self.errorLog('Problem with %s.tarZipSaveSpace()' % self.prog())
#########################################################################################################################
    ### <10> ### Extra Functions                                                                                        #
#########################################################################################################################
    def motifSeq(self):     ### Outputs sequence files for given motifs
        '''Outputs sequence files for given motifs.'''
        try:
            ### Setup ###
            if not self.dict['MotifSeq']: return False
            filedict = {}  # Motif:FileName
            self.obj['MotifSeq'] = rje_motiflist.MotifList(self.log,self.cmd_list+['motifs=None','minic=0','minpep=0','minfix=0'])
            for pattern in self.dict['MotifSeq']:
                Motif = self.obj['MotifSeq']._addMotif(pattern,pattern,reverse=False,check=False,logrem=True)
                if Motif and Motif.slimCode(): filedict[Motif] = self.dict['MotifSeq'][pattern]
                else: self.log.errorLog('Cannot make Motif object and/or slimCode for "%s"' % pattern,printerror=False)
            if not filedict: return True

            ### Sequences ###
            for Motif in filedict:
                ## Make SeqList ##
                mseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                slim = Motif.info['Slim']
                if Motif.info['Slim'] not in self.dict['Slim']:
                    try:
                        for Seq in self.seqs():
                            if Motif.searchSequence(Seq.info['MaskSeq']): mseq.seq.append(Seq)
                    except:
                        self.log.errorLog('SLiM Error (%s)' % slim,quitchoice=False)
                        continue
                else:
                    for (Seq,pos) in self.dict['Slim'][slim]['Occ']:
                        if Seq not in mseq.seq: mseq.seq.append(Seq)
                ## Save ##
                mseq.saveFasta(seqfile=filedict[Motif])
                del mseq
            return True
        except: self.log.errorLog('Major problem with %s.motifSeq()' % self.prog())
        return False
#########################################################################################################################
    def randomise(self):    ### Makes random datasets using batch file UPCs
        '''Makes random datasets using batch file UPCs.'''
        try:### ~ [1] ~ Setup Random sequence source ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFile(self.info['RandSource']):
                randsource = rje_seq.SeqList(self.log,['accnr=F','seqnr=F']+self.cmd_list+['seqin=%s' % self.info['RandSource'],'autoload=T'])
                if randsource.seqNum() < 1: return self.errorLog('Problem with %s' % self.info['RandSource'])
            else:
                randsource = None
                self.printLog('#RAND','No Random Sequence source - will use actual UPC')
            ### ~ [2] ~ Setup UP Lists and make *.upc if missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fullupc = []    # Full list of UPCs (containing Sequence objects!)
            datsize = []    # List of dataset sizes (UPC number) to make
            if self.list['Batch']:
                batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
            if not self.list['Batch'] or not batchfiles:
                self.log.errorLog('No input file(s) found!',printerror=False)
                return False
            ## ~ [2a] ~ Read *.upc and check/create *.slimdb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for infile in batchfiles:
                ## Basefile ##
                self.info['Basefile'] = rje.baseFile(infile)
                if self.info['ResDir'].lower() not in ['','none']:
                    rje.mkDir(self,self.info['ResDir'])
                    self.info['Basefile'] = self.info['ResDir'] + self.dataset()
                ## SeqNum ##
                seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['seqin=%s' % infile,'autoload=T']
                self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
                if self.seqNum() > self.stat['MaxSeq'] > 0: 
                    self.log.printLog('#SKIP','Skipping %s: %s seq > %d MaxSeq' % (self.info['Basefile'],self.seqNum(),self.stat['MaxSeq']))
                    continue
                ## RandSource/UPC ##
                if randsource:
                    datsize.append(self.seqNum())
                    for i in range(self.seqNum()):
                        r = random.randint(0,randsource.seqNum()-1)
                        fullupc.append((randsource.seq[r],))                        
                else:
                    self.makeUPC()
                    datsize.append(self.UPNum())
                    for upc in self.list['UP']: fullupc.append(upc)
            ### ~ [3] ~ Make Random Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fullupc = rje.randomList(fullupc)
            rje.mkDir(self,self.info['RanDir'])
            rbase = self.info['RanDir'] + self.info['Randbase']
            r = 0
            for d in datsize:
                r += 1
                dseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                for u in range(d):
                    upc = fullupc.pop(0)
                    for seq in upc: dseq.seq.append(seq)
                dseq._checkForDup(remdup=True)
                dseq.info['Name'] = '%s_%s_%d-%d.fas' % (rbase,rje.preZero(r,len(datsize)),dseq.seqNum(),d)
                dseq.saveFasta()
            ### ~ [4] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Win32'] and len(sys.argv): self.verbose(0,0,'Finished!',1)
        except:
            self.log.errorLog('Major problem with %s.randomise()' % self.prog())
            raise
#########################################################################################################################            
### End of SECTION II: SLiMCore Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def patternFromCode(slim):  ### Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG)
    '''Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG).'''
    pattern = ''
    slist = string.split(slim,'-')
    i = 0
    while i < (len(slist)):
        a = slist[i]
        i += 2
        x = '0'
        if i < len(slist): x = slist[i-1]
        if len(a) == 1: pattern += a 
        else: pattern += '[%s]' % a
        if len(x) == 1: pattern += '.' * int(x)
        else: pattern += '.{%d,%d}' % (int(x[0]),int(x[-1]))
    return pattern
#########################################################################################################################
def slimPos(slim):  ### Returns the number of positions in a slim
    return (string.count(slim,'-') / 2) + 1
#########################################################################################################################
def slimLen(slim):  ### Returns length of slim
    return len(patternFromCode(slim))
#########################################################################################################################
def slimDif(slim1,slim2):   ### Returns number of different positins between slim1 and slim2
    '''Returns number of different positions between slim1 and slim2.'''
    if slim1 == slim2: return 0
    (split1,split2) = (string.split(slim1,'-'),string.split(slim2,'-'))
    if len(split1) != len(split2): return rje.modulus(len(split1) - len(split2))
    d = 0
    for i in range(len(split1)):
        if split1[i] != split2[i]: d += 1
    return d
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try: SLiMCore(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: pass    #!#return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
