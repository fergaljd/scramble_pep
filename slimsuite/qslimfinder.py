#!/usr/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       qslimfinder
Description:  Query Short Linear Motif Finder
Version:      1.2
Last Edit:    15/09/10
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    QSLiMFinder is a modification of the basic SLiMFinder tool to specifically look for SLiMs shared by a query sequence
    and one or more additional sequences. To do this, SLiMBuild first identifies all motifs that are present in the query
    sequences before removing it (and its UPC) from the dataset. The rest of the search and stats takes place using the
    remainder of the dataset but only using motifs found in the query. The final correction for multiple testing is made
    using a motif space defined by the original query sequence, rather than the full potential motif space used by the
    original SLiMFinder. This is offset against the increased probability of the observed motif support values due to the
    reduction of support that results from removing the query sequence but could potentially still identify SLiMs will
    increased significance.

    Note that minocc and ambocc values *include* the query sequence, e.g. minocc=2 specifies the query and ONE other UPC.    
    
Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options ### 
    seqin=FILE      : Sequence file to search [None]
    batch=LIST      : List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    query=LIST      : Return only SLiMs that occur in 1+ Query sequences (Name/AccNum/Seq Number) [1]
    addquery=FILE   : Adds query sequence(s) to batch jobs from FILE [None]
    maxseq=X        : Maximum number of sequences to process [500]
    maxupc=X        : Maximum UPC size of dataset to process [0]
    sizesort=X      : Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
    walltime=X      : Time in hours before program will abort search and exit [1.0]
    resfile=FILE    : Main SLiMFinder results table [slimfinder.csv]
    resdir=PATH     : Redirect individual output files to specified directory (and look for intermediates) [SLiMFinder/]
    buildpath=PATH  : Alternative path to look for existing intermediate files [SLiMFinder/]
    force=T/F       : Force re-running of BLAST, UPC generation and SLiMBuild [False]
    pickup=T/F      : Pick-up from aborted batch run by identifying datasets in resfile using RunID [False]
    dna=T/F         : Whether the sequences files are DNA rather than protein [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SLiMBuild Options I: Evolutionary Filtering  ###
    efilter=T/F     : Whether to use evolutionary filter [True]
    blastf=T/F      : Use BLAST Complexity filter when determining relationships [True]
    blaste=X        : BLAST e-value threshold for determining relationships [1e=4]
    altdis=FILE     : Alternative all by all distance matrix for relationships [None]
    gablamdis=FILE  : Alternative GABLAM results file [None] (!!!Experimental feature!!!)
    homcut=X        : Max number of homologues to allow (to reduce large multi-domain families) [0]

    ### SLiMBuild Options II: Input Masking ###
    masking=T/F     : Master control switch to turn off all masking if False [True]
    dismask=T/F     : Whether to mask ordered regions (see rje_disorder for options) [False]
    consmask=T/F    : Whether to use relative conservation masking [False]
    ftmask=LIST     : UniProt features to mask out [EM]
    imask=LIST      : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
    compmask=X,Y    : Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    casemask=X      : Mask Upper or Lower case [None]
    motifmask=X     : List (or file) of motifs to mask from input sequences []
    metmask=T/F     : Masks the N-terminal M (can be useful if termini=T) [True]
    posmask=LIST    : Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    aamask=LIST     : Masks list of AAs from all sequences (reduces alphabet) []
    qregion=X,Y     : Mask all but the region of the query from (and including) residue X to residue Y [0,-1]
    
    ### SLiMBuild Options III: Basic Motif Construction ###
    termini=T/F     : Whether to add termini characters (^ & $) to search sequences [True]
    minwild=X       : Minimum number of consecutive wildcard positions to allow [0]
    maxwild=X       : Maximum number of consecutive wildcard positions to allow [2]
    slimlen=X       : Maximum length of SLiMs to return (no. non-wildcard positions) [5]
    minocc=X        : Minimum number of unrelated occurrences for returned SLiMs. (Proportion of UP if < 1) [0.05]
    absmin=X        : Used if minocc<1 to define absolute min. UP occ [3]
    alphahelix=T/F  : Special i, i+3/4, i+7 motif discovery [False]

    ### SLiMBuild Options IV: Ambiguity ###
    preamb=T/F      : Whether to search for ambiguous motifs during motif discovery [True]
    ambocc=X        : Min. UP occurrence for subvariants of ambiguous motifs (minocc if 0 or > minocc) [0.05]
    absminamb=X     : Used if ambocc<1 to define absolute min. UP occ [2]
    equiv=LIST      : List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMVF,FYW,FYH,KRH,DE,ST]
    wildvar=T/F     : Whether to allow variable length wildcards [True]
    combamb=T/F     : Whether to search for combined amino acid degeneracy and variable wildcards [False]

    ### SLiMBuild Options V: Advanced Motif Filtering ###
    musthave=LIST   : Returned motifs must contain one or more of the AAs in LIST (reduces search space) []
    focus=FILE      : FILE containing focal groups for SLiM return (see Manual for details) [None]
    focusocc=X      : Motif must appear in X+ focus groups (0 = all) [0]
    * See also rje_slimcalc options for occurrence-based calculations and filtering *
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SLiMChance Options ###
    slimchance=T/F  : Execute main SLiMFinder probability method and outputs [True]
    probcut=X       : Probability cut-off for returned motifs [0.1]
    maskfreq=T/F    : Whether to use masked AA Frequencies (True), or (False) mask after frequency calculations [False]
    aafreq=FILE     : Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    aadimerfreq=FILE: Use empirical dimer frequencies from FILE (fasta or *.aadimer.tdt) (!!!Experimental!!!) [None]
    negatives=FILE  : Multiply raw probabilities by under-representation in FILE (!!!Experimental!!!) [None]
    smearfreq=T/F   : Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    seqocc=T/F      : Whether to upweight for multiple occurrences in same sequence (heuristic) [False]
    probscore=X     : Score to be used for probability cut-off and ranking (Prob/Sig) [Sig]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Advanced Output Options I: Output data ###
    clouds=X        : Identifies motif "clouds" which overlap at 2+ positions in X+ sequences (0=minocc / -1=off) [2]
    runid=X         : Run ID for resfile (allows multiple runs on same data) [DATE:TIME]
    logmask=T/F     : Whether to log the masking of individual sequences [True]
    slimcheck=FILE  : Motif file/list to add to resfile output [] 

    ### Advanced Output Options II: Output formats ###
    teiresias=T/F   : Replace TEIRESIAS, making *.out and *.mask.fasta files [False]
    slimdisc=T/F    : Emulate SLiMDisc output format (*.rank & *.dat.rank + TEIRESIAS *.out & *.fasta) [False]
    extras=X        : Whether to generate additional output files (alignments etc.) [1]
                        - 0 = No output beyond main results file
                        - 1 = Generate occurrence file, alignments and cloud file
                        - 2 = Generate all additional SLiMFinder outputs
                        - 3 = Generate SLiMDisc emulation too (equiv extras=2 slimdisc=T)
    targz=T/F       : Whether to tar and zip dataset result files (UNIX only) [False]
    savespace=0     : Delete "unneccessary" files following run (best used with targz): [0]
                        - 0 = Delete no files
                        - 1 = Delete all bar *.upc and *.pickle files
                        - 2 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)

    ### Advanced Output Options III: Additional Motif Filtering ### 
    topranks=X      : Will only output top X motifs meeting probcut [1000]
    minic=X         : Minimum information content for returned motifs [2.1]
    * See also rje_slimcalc options for occurrence-based calculations and filtering *
    
Uses general modules: copy, glob, math, os, string, sys, time
Uses RJE modules: slimfinder, rje, rje_blast, rje_slim, rje_slimlist, rje_slimcalc, rje_slimcore, rje_dismatrix_V2,
    rje_seq, rje_scoring
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, math, os, pickle, string, sys, time
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_blast, rje_seq, rje_sequence, rje_scoring, rje_xgmml
import slimfinder, rje_slim, rje_slimcalc, rje_slimcore, rje_slimlist
import rje_motif_V3 as rje_motif            # Used for expect method only 
import rje_dismatrix_V2 as rje_dismatrix
import comparimotif_V3 as comparimotif
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on SLiMFinder 3.5.
    # 1.0 - Test & Modified to include AA masking.
    # 1.1 - Added sizesort.
    # 1.2 - Added the addquery function.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Finish and test.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('QSLiMFinder', '1.2', 'September 2010', '2008')
    description = 'Query Short Linear Motif Finder'
    author = 'Richard J. Edwards, Norman E. Davey & Denis C. Shields'
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
            if rje.yesNo('Show SLiMCalc commandline options?'): out.verbose(-1,4,text=rje_slimcalc.__doc__)
            if rje.yesNo('Show RJE_SEQ commandline options?'): out.verbose(-1,4,text=rje_seq.__doc__)
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
default_equiv = 'AGS,ILMVF,FYW,FYH,KRH,DE,ST'
basic_headers = ['Rank','Pattern','IC','Occ','Support','UP','ExpUP','Prob','Sig','Cloud','CloudSeq','CloudUP']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMFinder Class                                                                                        #
#########################################################################################################################
class QSLiMFinder(slimfinder.SLiMFinder):     
    '''
    QSLiMFinder Class. Author: Rich Edwards (2008).

    See SLiMFinder Class for details of Attributes.    
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    # def _setAttributes(self):   #!# Same defaults etc as SLiMFinder
    # def _cmdList(self): #!# Same commands as SLiMFinder
#########################################################################################################################
    ### <2> ### Simple Stats Methods                                                                                    #
#########################################################################################################################
    def QUPNum(self): return self.UPNum() - len(self.dict['FocusUPC']['Query'])
    def QUP(self,upc): return upc in self.dict['FocusUPC']['Query']
    def slimQUP(self,slim):     ### Returns number of UPC for a SLiM, excluding Query UPC
        if not self.dict['Slim'].has_key(slim): return 0
        qupx = 0
        for upc in self.dict['Slim'][slim]['UP']:
            if not self.QUP(upc): qupx += 1
        return qupx
#########################################################################################################################
    ### <3> ### General Run Methods                                                                                     #
#########################################################################################################################
    def run(self,batch=False):  ### Main SLiMFinder Run Method
        '''
        Main SLiMFinder Run Method:
        0. PreCheck:
            - Check for randomise function and execute if appropriate
        1. Input:
            - Read sequences into SeqList
            - or - Identify appropriate Batch datasets and rerun each with batch=True
        2. SLiMBuild:
            - Check for existing Pickle and load if found. Check appropriate parameter settings and re-run if desired.
            - or - Save sequences as fasta and mask sequences in SeqList
            -  Perform BLAST and generate UPC based on saved fasta.
            - Calculate AAFreq for each sequence and UPC.
            - Find all dimer motifs in dataset using MinWild/MaxWild parameters.
            - Extend to SLiMs and add ambiguity
        5. Identify significant SLiMs.
        6. Output results and tidy files.
        >> batch:bool [False] = whether this run is already an individual batch mode run.
        '''
        try:###~PRECHECK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ### Randomise Function ###
            if self.opt['Randomise'] and not batch: return self.randomise()

            ###~INPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            seqcmd = ['gnspacc=T','usecase=T','query=1'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            self.setupBasefile()
            self.loadAADimerFreq()
            ## Batch Mode ##
            if not batch and self.seqNum() < 1:   # No sequences loaded - use batch mode
                pickup = None
                if self.opt['Pickup'] and os.path.exists(self.info['ResFile']):
                    try:
                        #pickup = rje.readDelimit(self.loadFromFile(self.info['ResFile'])[-1],rje.delimitFromExt(filename=self.info['ResFile']))[0]
                        pickup = rje.dataDict(self,self.info['ResFile'],['RunID'],['Dataset'],lists=True)[self.info['RunID']]['Dataset']
                    except: self.log.errorLog('Could not read pickup dataset from %s' % self.info['ResFile'])
                #batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
                #self.log.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
                batchfiles = self.batchFiles(pickup)
                if not batchfiles: self.log.errorLog('No input files found!',printerror=False)
                else:
                    self.setupResults()                 ## Sets up OccStats filter etc. - check against Pickle ##
                    self.backupOrCreateResFile()
                    self.opt['Append'] = True
                    mycmd = self.cmd_list[0:]
                    self.list['Batch'] = []
                    bx = 0
                    for infile in batchfiles:
                        bx += 1
                        if pickup:
                            next = os.path.split(rje.baseFile(infile))[1]
                            #if pickup == os.path.split(rje.baseFile(infile))[1]: pickup = None
                            if next in pickup:
                                self.printLog('#PICKUP','Skipping %s batch file %s %s' % (self.info['RunID'],rje.integerString(bx),infile),log=False)
                                continue
                        self.printLog('#BATCH','Batch running %s' % infile)
                        bsf = self.newBatchRun(infile)
                        bsf.dict['AADimerFreq'] = self.dict['AADimerFreq']
                        bsf.run(batch=True)
                        self.log.printLog('#BATCH','Batch file %s run. Cleaning up for next file.' % infile)
                        del bsf.obj
                        del bsf.list
                        del bsf.dict
                        del bsf
                        self.log.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(bx),rje.integerString(len(batchfiles)-bx)),log=False)
                if self.opt['Win32'] and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return
            ## Check whether to bother running dataset at all - Check Input versus Min and Max Seq ##
            if self.stat['MaxSeq'] > 0 and self.stat['MaxSeq'] < self.obj['SeqList'].seqNum():
                self.log.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.integerString(self.obj['SeqList'].seqNum()),rje.integerString(self.stat['MaxSeq'])))
                return False
            if self.seqNum() < self.stat['MinOcc']:
                self.log.printLog('#SEQ','Insufficient Sequences (%d) for MinOcc setting (%d). Run aborted.' % (self.seqNum(),self.stat['MinOcc']))
                return

            ###~SLiMBuild~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.stat['StartTime'] = time.time()
            self.addQuery()
            ## UPC and MinOcc settings: needed to identify correct pickle so must be done first ##
            if not self.makeUPC():
                return self.log.errorLog('Error during makeUPC(). Abandoning %s run' % self.dataset(),printerror=False)
            if not self.setupFocus(): return    ## Setup Focus after UPC & MST - also check against Pickle ##
            if not self.setupQueryFocus(): return
            self.reportQueryUPC()
            if not self.setupMinOcc(): 
                self.printLog('#UPC','Insufficient UPC (%d) for MinOcc setting (%d). Run aborted.' % (self.UPNum(),self.stat['MinOcc']))
                return
            if self.stat['MaxUPC'] >= self.stat['MinOcc'] and self.stat['MaxUPC'] < self.UPNum():
                self.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). Run aborted.' % (self.UPNum(),self.stat['MaxUPC']))
                return

            ## Check for existing pickle to replace SLiMBuild portion ##
            self.setupResults()                 ## Sets up OccStats filter etc. - check against Pickle ##
            pickled = self.pickleMe(load=self.opt['Pickle'])  # Returns appropriate pickled SLiMFinder Object, else None
            if pickled: self = pickled  ## Replace me with my pickle!

            ## Setup Main Results File early in case of user intervention ##
            self.backupOrCreateResFile()

            ## AA Frequency Calculations made early as needed superficially in SLiMBuild ##
            if not pickled: self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
            if self.opt['MaskFreq']: self.makeAAFreq()
            else: 
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
                self.makeAAFreq()
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
            self.adjustAATotals()

            ## Execute SLiMBuild if pickle not loaded, else recalculate Bonferroni ##
            if self.opt['SlimBuild'] or self.opt['SlimChance'] or self.opt['SlimDisc']:
                if not pickled:     ## Find all dimer motifs in dataset using MaxWild parameters. ##
                    self.makeBonferroni()   # Estimates and reports total no. motifs in dataset
                    self.makeDimers()       # Makes all ai.{0,x}aj dimers
                    self.reduceDimers()     # Reduces to interesting subset
                    self.makeSLiMs()        # Makes all SLiMs with sufficient support
                    self.pickleMe()         # Generates pickling for speedy re-running
                else:
                    for i in range(2,self.stat['SlimLen']+1): self.printLog('#SPACE','Motif Space, %d positions: %s motifs' % (i,rje.integerString(self.dict['Extremf.'][i])))

            ### Special MotifSeq Output ###
            # This must occur after Input masking but needs no AA Frequencies or SLiMBuild #
            if self.motifSeq() and not self.opt['SlimBuild']: return     

            ###~Post-SLiMBuild Processing/Filtering before SLiMChance and Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
            ## Non-SLiMChance filtering of motifs ##
            self.dict['ElementIC'] = {}
            self.filterSLiMs()
            ## TEIRESIAS Output ##
            if self.opt['Teiresias'] or self.opt['SlimDisc']: self.teiresias()
            if not self.opt['SlimChance'] and not self.opt['SlimDisc']:
                self.log.printLog('#PROB','SlimChance=F and SlimDisc=F : No SLiM probability calculations')
                return

            ###~SLiMChance Probability and Significance Calculations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ptxt = 'Calculating SLiM Probabilities (p<=%s)' % rje_slim.expectString(self.stat['ProbCut'])
            sx = 0
            self.dict['PCut'] = {}
            for slim in self.dict['Slim']:
                self.log.printLog('\r#PROB','%s %.2f%% (%s Sig SLiMs)' % (ptxt,sx/self.slimNum(),rje.integerString(len(self.list['SigSlim']))),newline=False,log=False)
                sx += 100.0
                self.sigSlim(slim)
                self.wallTime()
            self.log.printLog('\r#PROB','%s complete: %s Sig SLiMs.' % (ptxt,rje.integerString(len(self.list['SigSlim']))))
            ## Reduce to Significant SLiMs ##
            self.obj['SlimList'].list['Motif'] = []
            for slim in self.dict['Slim'].keys()[0:]:
                if slim not in self.list['SigSlim']: self.dict['Slim'].pop(slim)
            for slim in self.list['SigSlim']: self.obj['SlimList'].list['Motif'].append(self.addSLiMToList(slim))

            ###~SLiMFinder Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            if self.list['SigSlim']: self.progLog('\r#OCC','Calculating Occ Stats etc...')
            ### Output results and tidy files. ###
            if not (self.occFilter() or self.statFilter()): self.calculateSLiMOccStats()
            self.obj['SlimList'].combMotifOccStats()
            self.tidyMotifObjects()     # Temporary solution to problem with unknown cause
            self.makeClouds()           # Identifies "clouds" of motifs - similar patterns that overlap            
            self.rankScore()            # Converts rankings into Numeric
            self.results()              # Controls SLiMFinder results output
            self.slimCheck()            # Additional SlimCheck Motifs 
            if self.extras(): self.extraOutput()   # MotifList Outputs 
            self.tarZipSaveSpace()      # Tarring, Zipping and Saving Space 

            ### End ###
            self.log.printLog('#RES','SLiMFinder results output to %s and %s*' % (self.info['ResFile'],self.info['Basefile']))
            if self.opt['Win32'] and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            return True
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.stat['WallTime'] <= 0 or (time.time() - self.stat['StartTime']) < (self.stat['WallTime']*3600): raise
            if self.list['Headers']: self.results(aborted=True)
            return False # Walltime reached
        except:
            self.log.errorLog('Error in SLiMFinder.run()',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def newBatchRun(self,infile):   ### Returns SLiMFinder object for new batch run
        '''Returns QSLiMFinder object for new batch run.'''
        return QSLiMFinder(self.log,self.cmd_list[0:] + ['seqin=%s' % infile,'append=%s' % self.opt['Append']])
#########################################################################################################################
    ### <4> ### Setup/Input Methods                                                                                     #
#########################################################################################################################
    def reportQueryUPC(self):   ### Reports input with UPC similarity to Query - also of interest
        '''Reports input with UPC similarity to Query - also of interest.'''
        for upc in self.dict['FocusUPC']['Query']:
            for seq in upc:
                if seq not in self.dict['Focus']['Query']:
                    self.printLog('#QUPC','%s is found in same UPC as Query' % seq.shortName())
#########################################################################################################################
    def addQuery(self):     ### Loads and sets addQuery sequence(s)
        '''Loads and sets addQuery sequence(s).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try:
                if self.info['AddQuery'].lower() in ['','none']: return
            except: return
            if not os.path.exists(self.info['AddQuery']):
                self.errorLog('AddQuery file "%s" not found!' % self.info['AddQuery'],printerror=False)
                raise IOError
            ### ~ [1] ~ Load and set queries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            qseq = rje_seq.SeqList(self.log,seqcmd+['seqin=%s' % self.info['AddQuery']])
            self.obj['SeqList'].seq += qseq.seq
            self.printLog('#ADD','%d sequence(s) added from %s' % (qseq.seqNum(),self.info['AddQuery']))
            if not self.list['Query']:
                self.list['Query'] = qseq.accList()
                self.printLog('#QRY','Query set: %s' % string.join(self.list['Query'],'; '))
        except: self.errorLog('QSLiMFinder.addQuery() error.'); raise
#########################################################################################################################
    ### <5> ### SLiMBuild Generation Methods                                                                            #
#########################################################################################################################
    def makeDimers(self):   ### Finds all possible dimers with wildcards, using MaxWild stat
        '''Finds all possible dimers with wildcards, using MinWild/MaxWild stat.'''
        try:
            ### Setup ###
            self.dict['Dimers'] = {}
            self.dict['DimFreq'] = {}
            uplist = self.list['UP'][0:]
            nonx = {}   # Total count of non-X positions in UPC
            for upc in uplist + self.seqs():
                nonx[upc] = 0.0
                self.dict['DimFreq'][upc] = [0] * (self.stat['MaxWild'] + 1)
            dx = 0
            sx = 0

            ### Read ###
            for seq in self.seqs():
                ## Setup UPC ##
                upc = self.getUP(seq)
                if upc in self.dict['FocusUPC']['Query'] and seq not in self.dict['Focus']['Query']: continue
                ## Setup Sequence ##
                sx += 1
                sequence = seq.info['Sequence'].upper()
                if self.opt['DNA']: sequence = string.replace(sequence,'N','X')
                if self.opt['Termini']: sequence = '^%s$' % sequence
                ## Find Dimers ##
                for i in range(len(sequence)):
                    ## Choose first position and check for wildcard ##
                    r = i
                    if self.opt['Termini']: r = i - 1
                    ai = sequence[i]
                    if ai in wildcards: continue
                    nonx[upc] += 1
                    nonx[seq] += 1
                    ## Examine each wildcard length in turn ##
                    for x in range(self.stat['MinWild'],self.stat['MaxWild']+1):
                        j = i + x + 1
                        if len(sequence) <= j: continue
                        aj = sequence[j]
                        if aj in wildcards: continue
                        ## Add Dimer ##
                        self.dict['DimFreq'][seq][x] += 1
                        self.dict['DimFreq'][upc][x] += 1
                        if not self.dict['Dimers'].has_key(ai): self.dict['Dimers'][ai] = {}
                        if not self.dict['Dimers'][ai].has_key(x): self.dict['Dimers'][ai][x] = {}
                        if not self.dict['Dimers'][ai][x].has_key(aj):
                            self.dict['Dimers'][ai][x][aj] = {'UP':[],'Occ':[]}
                            dx += 1
                        if upc not in self.dict['Dimers'][ai][x][aj]['UP']: self.dict['Dimers'][ai][x][aj]['UP'].append(upc)
                        self.dict['Dimers'][ai][x][aj]['Occ'].append((seq,r))
                        newslim = '%s-%s-%s' % (ai,x,aj)   #!# Str->List mod 1.5 #!#
                        if not self.dict['SeqOcc'].has_key(newslim): self.dict['SeqOcc'][newslim] = {seq:1}
                        elif not self.dict['SeqOcc'][newslim].has_key(seq): self.dict['SeqOcc'][newslim][seq] = 1
                        else: self.dict['SeqOcc'][newslim][seq] += 1
                self.log.printLog('\r#DIM','Reading dimers (%d seq) %s dimers' % (sx,rje.integerString(dx)),newline=False,log=False)
            self.log.printLog('\r#DIM','Read dimers from %d seq: %s dimers' % (sx,rje.integerString(dx)))
            self.stat['Dimers'] = dx

            ### Adjust DimFreq ###
            for upc in uplist + self.seqs():
                for x in range(self.stat['MinWild'],(self.stat['MaxWild'] + 1)):
                    if nonx[upc]: self.dict['DimFreq'][upc][x] = self.dict['DimFreq'][upc][x] / nonx[upc]
                    else:
                        if upc in self.list['UP']: self.printLog('#NONX','WARNING! UPC cluster %d has zero unmasked residues!' % uplist.index(upc))
                        elif self.getUP(upc) not in self.dict['FocusUPC']['Query'] or upc in self.dict['Focus']['Query']: self.printLog('#NONX','WARNING! Sequence %s has zero unmasked residues!' % upc.shortName())
                        self.dict['DimFreq'][upc][x] = 0.0

        except:
            self.log.errorLog('Major problem during makeDimers()')
            raise
#########################################################################################################################
    def reduceDimers(self):     ### Reduces Dimers to those with enough Support
        '''Reduces Dimers to those with enough Support.'''
        try:### ~ [1] ~ Select Dimers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dx = 0      # Dimer count with sufficient Occ
            qdx = 0     # Total dimer count in Query
            for ai in self.dict['Dimers'].keys()[0:]:
                for x in self.dict['Dimers'][ai].keys()[0:]:
                    for aj in self.dict['Dimers'][ai][x].keys()[0:]:
                        inquery = False
                        for upc in self.dict['Dimers'][ai][x][aj]['UP']: inquery = inquery or upc in self.dict['FocusUPC']['Query']
                        if inquery: qdx += 1
                        ox = len(self.dict['Dimers'][ai][x][aj]['UP'])
                        # Reject dimer not in query and no good for ambiguity
                        if ox < self.stat['AmbOcc'] and not inquery: self.dict['Dimers'][ai][x].pop(aj)
                        else: dx += 1
                    self.progLog('\r#DIM','Reducing dimers: %s >= %d of %d UPC ' % (rje.integerString(dx),self.stat['AmbOcc'],self.UPNum()))
                if not self.dict['Dimers'][ai][x]: self.dict['Dimers'][ai].pop(x)
            if not self.dict['Dimers'][ai]: self.dict['Dimers'].pop(ai)
            self.printLog('\r#DIM','Reducing dimers: %s >= %d of %d UPC ' % (rje.integerString(dx),self.stat['AmbOcc'],self.UPNum()))
            ### ~ [2] ~ Reduce Bonferroni ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SPACE','Query Motif Space, 2 positions: %s motifs -> %s query motifs' % (rje.integerString(self.dict['Extremf.'][2]),rje.integerString(qdx)))
            self.dict['Extremf.'][2] = qdx
        except: self.errorLog('Problem reducing Dimers to AmbOcc+')
#########################################################################################################################
    def makeSLiMs(self):    ### Makes SLiMs with enough support from Dimers
        '''Makes SLiMs with enough support from Dimers.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Slim'] = {}
            prevslim = []
            ## ~ [0a] ~ Select Dimers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for ai in self.dict['Dimers']:
                for x in self.dict['Dimers'][ai]:
                    for aj in self.dict['Dimers'][ai][x].keys()[0:]:
                        slim = '%s-%s-%s' % (ai,x,aj)
                        prevslim.append(slim)
                        self.dict['Slim'][slim] = self.dict['Dimers'][ai][x][aj]
                    self.progLog('\r#SLIM','Selecting 2aa SLiMs: %s >= %d of %d UPC ' % (rje.integerString(self.slimNum()),self.stat['AmbOcc'],self.UPNum()))
            self.printLog('\r#SLIM','Selecting 2aa SLiMs: %s >= %d of %d UPC' % (rje.integerString(self.slimNum()),self.stat['AmbOcc'],self.UPNum()))
            ## ~ [0b] ~ Ambiguity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.ambSLiM(prevslim)
            ### ~ [1] ~ Extend ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for f in range(3,self.stat['SlimLen']+1):
                if not prevslim: break
                newslim = []
                ex = 0.0
                for slim in prevslim:
                    self.progLog('\r#SLIM','Extending %daa SLiMs >= %d of %d UPC: %.1f%%' % (f,self.stat['AmbOcc'],self.UPNum(),ex/len(prevslim)))
                    ex += 100.0
                    newslim += self.extendSLiM(slim)    # Why was there "if self.opt['PreAmb']" ???
                    self.wallTime()
                prevslim = newslim                
                self.printLog('\r#SLIM','Extending %daa SLiMs >= %d of %d UPC: %s SLiMs' % (f,self.stat['AmbOcc'],self.UPNum(),rje.integerString(len(prevslim))))
                ## Query Bonferroni ##
                qsx = 0
                for slim in newslim:
                    inquery = False
                    for upc in self.dict['Slim'][slim]['UP']: inquery = inquery or upc in self.dict['FocusUPC']['Query']
                    if inquery: qsx += 1
                self.printLog('#SPACE','Query Motif Space, %d positions: %s motifs -> %s query motifs' % (f,rje.integerString(self.dict['Extremf.'][f]),rje.integerString(qsx)))
                self.dict['Extremf.'][f] = qsx
                ## Add ambiguity ##                
                self.ambSLiM(prevslim)
            self.log.printLog('#SLIM','%s SLiMs >= %d of %d UPC' % (rje.integerString(self.slimNum()),self.stat['AmbOcc'],self.UPNum()))
        except:
            self.log.errorLog('Fatal error making SLiMs from Dimers')
            raise
#########################################################################################################################
    def extendSLiM(self,slim):  ### Finds and returns extensions of SLiM with sufficient support
        '''
        Finds and returns extensions of SLiM with sufficient support.
        >> slim:str = SLiM to extend (using dimers)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ai = slim[-1]
            if not self.dict['Dimers'].has_key(ai): return []
            extend = []
            slimocc = self.dict['Slim'][slim]['Occ']
            ### ~ [1] ~ Try dimers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for x in self.dict['Dimers'][ai]:
                if self.opt['AlphaHelix'] and slim.find('-%d-' % x) > 0: continue   # Not i,i+3/4,i+7
                for aj in self.dict['Dimers'][ai][x]:   # No longer check minocc as removed in prev makeSLiMs
                    newslim = slim + '-%s-%s' % (x,aj)   #!# Str->List mod 1.5 #!#
                    newocc = []
                    newup = []
                    inquery = False
                    for (seq,pos) in slimocc:
                        if (seq,pos+slimLen(slim)-1) in self.dict['Dimers'][ai][x][aj]['Occ']:
                            newocc.append((seq,pos))
                            upc = self.getUP(seq)
                            if upc not in newup: newup.append(upc)
                        inquery = inquery or seq in self.dict['Focus']['Query']
                    if len(newup) >= self.stat['AmbOcc'] or inquery:
                        extend.append(newslim)
                        self.dict['Slim'][newslim] = {'Occ':newocc,'UP':newup}
                        for (seq,pos) in newocc:
                            if not self.dict['SeqOcc'].has_key(newslim): self.dict['SeqOcc'][newslim] = {seq:1}
                            elif not self.dict['SeqOcc'][newslim].has_key(seq): self.dict['SeqOcc'][newslim][seq] = 1
                            else: self.dict['SeqOcc'][newslim][seq] += 1
            ### Return ###
            return extend
        except:
            try: self.log.errorLog('Problem extending SLiM "%s"' % slim,quitchoice=True)
            except: raise KeyboardInterrupt
        return []
#########################################################################################################################
    ### <6> ### SLiM Filtering Methods                                                                                  #
#########################################################################################################################
    def setupQueryFocus(self):  ### Sets up Query based on Focus dictionary - needed for QSLiMFinder run.
        '''
        Sets up Query based on Focus dictionary - needed for QSLiMFinder run.
        Returns True if OK, else False (which cancels run).
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Query' in self.dict['Focus']: return True   # Query present - all OK.
            ### ~ [1] ~ Take Query from other Focus group? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for grp in rje.sortKeys(self.dict['Focus']):
                if self.i() < 0 or rje.yesNo('Use "%s" Focus Group as query?' % grp):
                    self.dict['Focus']['Query'] = self.dict['Focus'].pop(grp)
                    self.dict['FocusUPC']['Query'] = self.dict['FocusUPC'].pop(grp)
                    return True
            while 'query=None' in self.obj['SeqList'].cmd_list: self.obj['SeqList'].cmd_list.remove('query=None')
            qry = self.obj['SeqList'].querySeq()
            if qry: self.dict['Focus']['Query'] = [qry]; return True
            self.errorLog('Query needed for QSLiMFinder. Quitting run.',printerror=False)
        except: self.errorLog('Problem setting up Query for %s' % self.dataset())
        return False
#########################################################################################################################
    def slimFocus(self,slim):   ### Returns True if slim if Focal sequence groups, else False
        '''Returns True if slim if Focal sequence groups, else False.'''
        ### ~ Query Focus ~ ###
        inquery = False
        for (seq,occ) in self.dict['Slim'][slim]['Occ']: inquery = inquery or seq in self.dict['Focus']['Query']
        if not inquery: return False
        elif self.stat['FocusOcc'] < 2: return True
        ### ~ Old SLiM Focus ~ ###
        maxfail = 0
        if self.stat['FocusOcc'] > 0: maxfail = len(self.dict['Focus']) - self.stat['FocusOcc']
        slimgrp = self.dict['Focus'].keys()     # Groups not accounted for
        for (seq,occ) in self.dict['Slim'][slim]['Occ']:
            for grp in slimgrp[0:]:
                if seq in self.dict['Focus'][grp]: slimgrp.remove(grp)
        if len(slimgrp) > maxfail: return False    # Too many group(s) not accounted for by occs
        return True
#########################################################################################################################
    ### <7> ### SLiMChance Probability Methods                                                                          #
#########################################################################################################################
    def slimProb(self,slim): ### Calculate Probabilities for given SLiM
        '''Calculate Probabilities for given SLiM.'''
        try:
            ###~Calculate prob of 1+ occ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            p1 = {}         # Dictionary of {upc:chance of 1+ occ in upc}
            ##~~Setup pattern and variable-lenght multiplier~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            poslist = []    # List of AA positions in SLiM
            wildlist = []   # List of wildcard lengths in SLiM
            wild = False    # Whether next part is a wildcard length
            mult = 1        # Variable-length multiplier
            minslimlen = 0  # Minimum SLiM length
            for part in string.split(slim,'-'):      # Split SLiM code in components
                ## Update lists ##
                if wild: wildlist.append(part)
                else: poslist.append(part); minslimlen += 1
                ## Calculate multiplier ##
                if wild:
                    (minx,maxx) = (self.stat['MaxWild'],0)
                    for x in part:
                        minx = min(minx,int(x))
                        maxx = max(maxx,int(x))
                    mult *= (int(maxx) - int(minx) + 1)
                    minslimlen += minx
                wild = not wild
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            for upc in self.list['UP']:
                if self.QUP(upc): continue  # Ignore Query UPC
                if self.dict['AADimerFreq']: (k,N,p) = self.aaDp1(slim,upc)
                else:
                    ## Setup  parameters for binomial ##
                    N = self.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
                    p = 1.0                                 # Probability of SLiM at each position
                    k = 1                                   # Number of successful trials (occurrences)
                    if self.opt['SeqOcc'] and self.slimOccNum(slim,upc) > 1: k = self.slimOccNum(slim,upc)
                    ## Calculate p and N from AAFreq and DimFreq ##
                    for pos in poslist:     # AA position
                        posfreq = 0.0
                        for aa in pos: posfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        p *= posfreq
                    if self.opt['DimFreq']:
                        for dim in wildlist:    # DimerFreq
                            dimfreq = 0.0
                            for x in dim:
                                try: dimfreq += self.dict['DimFreq'][upc][int(x)]   # Options for wildcard length
                                except: pass
                            N *= (dimfreq / len(dim))       # Mutliply by mean dimer frequency
                    else: N -= ((minslimlen-1) * self.dict['MST'][upc])
                    N = max(0,N)
                    N *= mult       # Each length variant is effectively another position the SLiM could occur
                    if p > 1: p = 1.0   # Cannot in reality have p > 1!
                    ## Calculate binomial ##
                    p1[upc] = rje.binomial(k,N,p,usepoisson=False,callobj=self)
            ## Extra verbosity. Remove at some point? ##
            self.verbose(2,3,'%s: %s' % (patternFromCode(slim),p1.values()),1)
            self.verbose(2,3,'%s: %s vs %s\n' % (patternFromCode(slim),self.slimUP(slim),sum(p1.values())),2)

            ###~Calculate overall probability of observed support~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## All observed occurrences ##
            self.dict['Slim'][slim]['ExpUP'] = sum(p1.values())    # Expected number of observed UPCs
            (k,n,p) = (self.slimQUP(slim), self.QUPNum(), self.dict['Slim'][slim]['ExpUP']/self.QUPNum()) # Use mean p1+
            if k <= 0: self.dict['Slim'][slim]['Prob'] = 1.0
            else: self.dict['Slim'][slim]['Prob'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
            ## Catch for binomial problems. Should no longer happen. ##
            if self.dict['Slim'][slim]['Prob'] <= 0:    # Shouldn't happen now! #
                self.log.errorLog('Probability for %s <= 0.0 due to numerical limitations: Given arbitrary 1e-16.!' % (patternFromCode(slim)),printerror=False)
                self.dict['Slim'][slim]['Prob'] = 1e-16
            ## Correction for restricted focal sequences ##
            self.focusAdjustment(slim)
            ###~Old Score calculations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.dict['Slim'][slim]['S'] = self.slimUP(slim) * self.slimIC(slim)
            self.dict['Slim'][slim]['R'] = self.dict['Slim'][slim]['S'] * self.slimUP(slim) / self.dict['Slim'][slim]['ExpUP']
        except:
            self.log.errorLog('Error with slimProb(%s)' % slim)
            self.dict['Slim'][slim]['Prob'] = 1.0
#########################################################################################################################
    def aaDp1(self,slim,upc):  ### Setup  parameters for p1+ binomial using AADimerFreq
        '''Setup  parameters for p1+ binomial.'''
        ### ~ [1] Setup Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        N = self.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
        p = 0.0                                 # Probability of SLiM at each position
        k = 1                                   # Number of successful trials (occurrences)
        if self.opt['SeqOcc'] and self.slimOccNum(slim,upc) > 1: k = self.slimOccNum(slim,upc)

        ### ~ [2] Special AADimerFreq option ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for slimvar in rje.listCombos(string.split(slim,'-')):
            v = rje.getFromDict(self.dict['AAFreq'][upc],slimvar[0],returnkey=False,default=0.0)
            slist = slimvar[0:]
            while slist:
                [i,x,j] = slist[:3]
                slist = slist[2:]
                try: v *= self.dict['AADimerFreq'][i][x][j]
                except: v = 0.0
            p += v
        return (k,N,p)
#########################################################################################################################
    def focusAdjustment(self,slim): ### Adjust raw probabilities according to focus dictionary
        '''
        Adjust raw probabilities according to focus dictionary.
        >> slim:str = SLiM for probability adjustment
        '''
        try:### ~ [1] Calculate probabilities for each focus group ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(self.dict['Focus']) > 1: return self.printLog('#FOCUS','Focus adjustment not yet implemented in QSLiMFinder!')
            if not self.dict['Focus']: return
            pgroup = {}
            for grp in self.dict['Focus']:
                a = len(self.dict['FocusUPC'][grp])     # No. of UPC in focal group
                b = self.slimUP(slim)                   # No. of UPC that the SLiM occurs in
                N = self.UPNum()                        # Total number of UPC
                pgroup[grp] = self.abNprob(a,b,N,1,'more')

            ### ~ [2] Adjust probability ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            meanp = rje.meanse(pgroup.values())[0]
            self.dict['Slim'][slim]['Prob'] *= rje.binomial(self.stat['FocusOcc'],len(pgroup),meanp,callobj=self)
        except:
            self.log.errorLog('Major problem with SLiMFinder.focusAdjustment()')
            raise
#########################################################################################################################
    ### <8> ### SLiMBuild Pickle Methods                                                                                #
#########################################################################################################################
    def myPickle(self):  ### Returns pickle identifier, also used for Outputs "Build" column (self.info['Build'])
        '''Returns pickle identifier, also used for Outputs "Build" column.'''
        ## Pickle name ##
        if self.info['Build'].lower() in ['','none']:
            amb = 0
            if self.opt['PreAmb'] and self.list['Equiv']: amb += 1
            if self.opt['PreAmb'] and self.opt['WildVar']: amb += 2
            #x#self.deBug('%s: %d' % (self.list['Equiv'],amb))
            if amb == 3 and self.opt['CombAmb']: amb = 4
            elif self.opt['AlphaHelix']: self.info['Build'] = 'alpha-o%da%d' % (self.stat['AmbOcc'],amb)
            elif self.stat['MinWild']: self.info['Build'] = 'l%dw%d-%do%da%d' % (self.stat['SlimLen'],self.stat['MinWild'],self.stat['MaxWild'],self.stat['AmbOcc'],amb)
            else: self.info['Build'] = 'l%dw%do%da%d' % (self.stat['SlimLen'],self.stat['MaxWild'],self.stat['AmbOcc'],amb)
        return '%s.%s' % (self.info['Build'],self.maskText())
#########################################################################################################################
    def pickleMe(self,load=False):  ### Loads existing pickle, or saves pickle for later!
        '''Saves pickle for later!.'''
        try:
            ### Setup ###
            if load and self.opt['Force']: return None      # Re-run SLiMBuild!
            elif not load and self.opt['MemSaver'] and not self.extras(): return None  # Do not save pickle
            ## Pickle name ##
            mypickle = self.myPickle()

            ### Load ###
            if load:
                newme = None    ## New SLiMFinder object, loaded from Pickle
                ## Check for file and load ##
                for pfile in [self.info['Basefile'],'%s%s' % (self.info['BuildPath'],self.dataset())]:
                    if not self.opt['Win32'] and os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)):
                        if os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                            if rje.isYounger('%s.%s.pickle.gz' % (pfile,mypickle),'%s.%s.pickle' % (pfile,mypickle)) == '%s.%s.pickle.gz' % (pfile,mypickle):
                                os.unlink('%s.%s.pickle' % (pfile,mypickle))
                        if not os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                            try: os.system('gunzip %s.%s.pickle.gz' % (pfile,mypickle))
                            except: self.log.errorLog('Cannot unzip %s.%s.pickle.gz' % (pfile,mypickle))
                    if os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                        self.log.printLog('#LOAD','Attempting to load SLiMFinder pickle.',log=False)
                        newme = pickle.load(open('%s.%s.pickle' % (pfile,mypickle),'r'))
                        self.log.printLog('#LOAD','SLiMFinder intermediate loaded: %s.%s.pickle.' % (pfile,mypickle))
                        if not self.opt['Win32']:
                            try:
                                if os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)): os.unlink('%s.%s.pickle.gz' % (pfile,mypickle))
                                os.system('gzip %s.%s.pickle' % (pfile,mypickle))
                                self.log.printLog('#GZIP','SLiMFinder %s.%s.pickle zipped.' % (pfile,mypickle))
                            except: self.log.errorLog('Cannot gzip %s.%s.pickle' % (pfile,mypickle))
                        break
                if not newme: return None
                ## Check other pertinent attributes - masking and additional filtering ##
                ## Note that MustHave and OccFilter filtering currently occur *after* SLiMBuild only ##
                changes = []
                for var in ['CompMask','CaseMask','MotifMask']:         # Info
                    if self.info[var] != newme.info[var]: changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['Masking','DisMask','ConsMask','MaskM']:   # Opt
                    if self.opt[var] != newme.opt[var]: changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['FTMask','IMask','Equiv']:      # List   
                    slist = self.list[var][0:]
                    nlist = newme.list[var][0:]
                    slist.sort()
                    nlist.sort()
                    if slist != nlist:
                        changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                if newme.list['Query'] != self.list['Query'] or self.info['Focus'] != self.info['Focus']:
                    self.printLog('#PICKLE','Query/Focus Parameters changed. Making new pickle.')
                    return None
                if newme.dict['Focus']:     # Post-Focus filtering is OK! Only worry if pickle has focus filtering
                    if rje.sortKeys(newme.dict['Focus']) != rje.sortKeys(self.dict['Focus']): # Assume all else is the same! #
                        changes.append(self.log.errorLog('Warning: "Focus" parameter mismatch', printerror=False, nextline=False))
                if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']:
                    if newme.list['QRegion'] and 'Focus' in newme.dict and 'Query' in newme.dict['Focus']:
                        if newme.list['QRegion'] != self.list['QRegion']:
                            changes.append(self.errorLog('Warning: "QRegion" parameter mismatch (%s vs %s)' % (newme.list['QRegion'],self.list['QRegion']), printerror=False, nextline=False))
                    else: changes.append(self.errorLog('Warning: "QRegion" parameter mismatch', printerror=False, nextline=False))
                ## Recreate or use pickle but add new commands ##
                #x#self.deBug(changes)
                if changes and (self.stat['Interactive'] < 0 or rje.yesNo('%d SLiMBuild parameter mismatches with pickle. Create new pickle?' % len(changes))):
                    self.log.printLog('#PICKLE','Parameters changed. Making new pickle.')
                    return None
                self.list['Warning'] += changes
                newme.cmd_list = self.cmd_list
                newme.setInfo(self.info)
                newme.setStat(self.stat)
                self.opt['Masked'] = newme.opt['Masked']
                self.opt['DNA'] = newme.opt['DNA']
                newme.setOpt(self.opt)
                newme.info['ResFile'] = self.info['ResFile']
                newme.info['ResDir'] = self.info['ResDir']
                newme.info['BuildPath'] = self.info['BuildPath']
                newme.stat['StartTime'] = self.stat['StartTime']
                newme.obj['SlimList'] = self.obj['SlimList']    # Should take SLiMCalc with it
                newme.setLog(self.log)
                self.setLog(self.log)    #x# newme.log.stat['StartTime'] = self.log.stat['StartTime']
                for mylist in ['Headers','MustHave','NewScore']: newme.list[mylist] = self.list[mylist]
                for mydict in ['Focus','NewScore']: newme.dict[mydict] = self.dict[mydict]
                newme.setupFocus()  #!# Need to convert to new Seq Objects - clean up at some point! #!#
                return newme
                
            ### Save ###
            if 'Pickle' in self.opt and not self.opt['Pickle']:
                self.log.printLog('#PICKLE','SLiMFinder pickling disabled with pickle=F.',log=True)
                return None
            self.log.printLog('#SAVE','Attempting to save SLiMFinder with pickle.',log=False)
            pickle.dump(self,open('%s.%s.pickle' % (self.info['Basefile'],mypickle),'w'))
            self.log.printLog('#SAVE','SLiMFinder intermediate saved as %s.%s.pickle (Python pickle).' % (self.info['Basefile'],mypickle))
            if not self.opt['Win32']:
                try:
                    pfile = self.info['Basefile']
                    if os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)): os.unlink('%s.%s.pickle.gz' % (pfile,mypickle))
                    os.system('gzip %s.%s.pickle' % (pfile,mypickle))
                    self.log.printLog('#GZIP','SLiMFinder %s.%s.pickle zipped.' % (pfile,mypickle))
                except: self.log.errorLog('Cannot gzip %s.%s.pickle' % (pfile,mypickle))
            return None
        
        except:
            self.log.errorLog('Major problem with SLiMFinder pickling!')
            return None
#########################################################################################################################
    ### <9> ### Results Output Methods                                                                                  #
#########################################################################################################################
    def getSlimProb(self,slim):   ### Returns appropriate SLiM Probability given settings
        '''Returns appropriate SLiM Score given settings.'''
        return self.dict['Slim'][slim][self.info['ProbScore']]
#########################################################################################################################
    def makeClouds(self):   ### Identifies "clouds" of motifs - similar patterns that overlap  (from SigSlim)
        '''Identifies "clouds" of motifs - similar patterns that overlap (from SigSlim).'''
        try:
            ###~Setup~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.dict['Clouds'] = {}
            ctxt = 'Making "Motif Clouds" for %d Sig Motifs' % len(self.list['SigSlim'])
            if self.stat['Clouds'] < 0: return self.log.printLog('#CLOUD','Not %s (clouds=X < 0)' % ctxt)
            if self.stat['Clouds'] < 1: self.stat['Clouds'] = self.stat['MinOcc']
            clouds = {}                             # Dictionary of slim:[slim cloud list] used in construction
            ctot = len(self.list['SigSlim']) * 3    # Attributes for Clouding progress
            lx = 0.0
            ##~Rework Occ in searchable dictionary with every defined position~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            occ = {}        # New dictionary of {slim:{seq:[pos]}}
            for slim in self.list['SigSlim']:
                self.log.printLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot),newline=False,log=False)
                lx += 100.0
                occ[slim] = {}
                for seq in self.seqs(): occ[slim][seq] = []     # List of defined positions in SLiM-seq pair
                for var in rje.listCombos(string.split(slim,'-')):
                    varslim = string.join(var,'-')
                    if varslim not in self.dict['AllSlims']: continue
                    for (seq,pos) in self.dict['AllSlims'][varslim]['Occ']:
                        r = pos     # Position offset
                        w = [0]     # List of variable wildcard offsets
                        slist = string.split(varslim,'-')
                        if r not in occ[slim][seq] and slist[0] != '^': occ[slim][seq].append(r)
                        i = 0   # Code index
                        while i < (len(slist)):
                            i += 2          # Next defined position
                            x = '0'         # Content of wildcard position
                            if i < len(slist): x = slist[i-1]
                            r += 1 + int(x[0])
                            ## Make new wildcard offsets ##
                            old = w[0:]
                            w = []
                            for b in old:
                                for n in x: w.append(b+int(n)-int(x[0]))
                            ## Update occ ##
                            if i >= len(slist) and slist[-1] == '$': continue
                            for n in w:
                                if (r+n) not in occ[slim][seq]: occ[slim][seq].append(r+n)
                #X#print slim, occ[slim]

            #!#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#!#
            #!# Note that, currently, variable wildcard lengths are not handled properly and the method doesn't try to#!#
            #!# work out WHICH wildcard length applies for any given occurrence, it just assumes any of them are fine.#!#
            #!# This will only be a problem in the rare occasion that two motifs physically overlap but NOT at their  #!#
            #!# defined positions and either or both of them have flexible wildcards. Rare enough, I hope!            #!#
            #!#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#!#
                
            ###~Make cloud partnerships~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            for slim1 in self.list['SigSlim']:
                self.log.printLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot),newline=False,log=False)
                lx += 100.0
                clouds[slim1] = []
                for slim2 in self.list['SigSlim']:
                    if slim1 == slim2: continue
                    #x#if slim2 in clouds: continue    # Already checked against rest
                    cx = 0  # Shared sequence count
                    for seq in occ[slim1]:
                        #X#print slim1, slim2, occ[slim1][seq],
                        #X#self.deBug(occ[slim2][seq])
                        if not occ[slim2][seq]: continue    # No occ for this seq
                        px = 0
                        for pos in occ[slim1][seq]:
                            if pos in occ[slim2][seq]: px += 1
                        if px >= 2: cx += 1
                        #X#print slim1, 'v', slim2, seq, px, cx, '>=', self.stat['Clouds']
                        if cx >= self.stat['Clouds']:    # We have a match!
                            clouds[slim1].append(slim2)
                            break   # Check no more sequences
            cbug = ''
            for slim1 in self.list['SigSlim']:
                cbug += patternFromCode(slim1) + ': '
                for slim2 in clouds[slim1]: cbug += patternFromCode(slim2) + '; '
                cbug += '\n'
            #self.deBug(cbug)

            ###~Make Actual Clouds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            cnum = 0
            for slim in self.list['SigSlim'][0:]:
                ## Make a cloud for each Sig SLiM not already in one ##
                self.log.printLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot),newline=False,log=False)
                lx += 100.0
                if slim not in clouds: continue     # Absorbed into another cloud
                ## Make a new cloud with this SLiM as the "best" variant ## 
                cnum += 1   # New cloud
                self.dict['Clouds'][cnum] = {'Best':patternFromCode(slim),'Seq':[],'Slim':[patternFromCode(slim)],'Sig':rje_slim.expectString(self.dict['Slim'][slim]['Sig'])}
                for seq in occ[slim]:
                    if occ[slim][seq]: self.dict['Clouds'][cnum]['Seq'].append(seq)
                self.dict['Slim'][slim]['Cloud'] = cnum
                cx = 0
                while cx != len(clouds[slim]):      # Keep absorbing and reducing clouds till no more match
                    cx = len(clouds[slim])
                    for slim2 in clouds[slim][0:]:
                        if slim2 == slim: continue
                        if slim2 not in clouds: continue    # Already taken
                        self.dict['Slim'][slim2]['Cloud'] = cnum
                        self.dict['Clouds'][cnum]['Slim'].append(patternFromCode(slim2))
                        for seq in occ[slim2]:
                            if occ[slim2][seq] and seq not in self.dict['Clouds'][cnum]['Seq']: self.dict['Clouds'][cnum]['Seq'].append(seq)
                        for newslim in clouds.pop(slim2):   # Combine clouds
                            if newslim not in clouds[slim]: clouds[slim].append(newslim)
                            #X# self.dict['Slim'][newslim]['Cloud'] = cnum
                ## Convert Cloud seqs into UPC ##
                self.dict['Clouds'][cnum]['UPC'] = []
                for seq in self.dict['Clouds'][cnum]['Seq']:
                    upc = self.getUP(seq)
                    if upc not in self.dict['Clouds'][cnum]['UPC']: self.dict['Clouds'][cnum]['UPC'].append(upc)
                    
            ### Finish ###
            #x#self.deBug(self.dict['Clouds'])
            self.log.printLog('\r#CLOUD','%s: %d Clouds' % (ctxt,cnum))
        except: self.log.errorLog('Problem %s' % ctxt)
#########################################################################################################################       
    def rankScore(self):  ### Scores and Ranks Sig Motifs
        '''Scores and Ranks Sig Motifs.'''
        try:
            ### Assign numerical rankings ###
            if not self.list['SigSlim']: return
            (prev,rank,prevp) = (1,1,1)
            self.log.printLog('\r#RANK','Rank calculations...',newline=False,log=False)
            for slim in self.list['SigSlim']:
                ## Assign Rank ##
                if self.opt['SlimDisc'] or self.dict['Slim'][slim][self.info['ProbScore']] > prev or self.dict['Slim'][slim]['Prob'] != prevp:
                    self.dict['Slim'][slim]['Rank'] = self.list['SigSlim'].index(slim) + 1
                else: self.dict['Slim'][slim]['Rank'] = rank
                (prev,rank,prevp) = (self.dict['Slim'][slim][self.info['ProbScore']],self.dict['Slim'][slim]['Rank'],self.dict['Slim'][slim]['Prob'])
            self.log.printLog('\r#RANK','Rank calculations complete')
        except:
            self.log.errorLog('Major disaster during SLiMFinder.rankScore()')
            raise
#########################################################################################################################
    def backupOrCreateResFile(self):    ### Backups up and/or creates main results file
        '''Backups up and/or creates main results file.'''
        if self.info['ResFile'].lower() not in ['','none'] and self.opt['SlimChance']:
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.info['ResFile'],write=True))
            rje.delimitedFileOutput(self,self.info['ResFile'],self.resHead(),delimit,rje_backup=True)
#########################################################################################################################
    def resHead(self):  ### Returns main Output headers
        '''Returns main Output headers.'''
        heads = ['Dataset','RunID','Masking','Build','Chance','RunTime','SeqNum','UPNum','AANum','MotNum'] + self.list['Headers']
        if self.opt['Test']: heads.insert(heads.index('Sig'),'E')
        return heads
#########################################################################################################################
    def getSigSlim(self,pattern):   ### Returns slimcode for given pattern (should be in SigSlim)
        '''Returns slimcode for given pattern (should be in SigSlim).'''
        for slim in self.list['SigSlim']:
            if patternFromCode(slim) == pattern: return slim
        return None
#########################################################################################################################
    def slimCheck(self):    ### Checks given list of Motifs and adds to output
        '''Checks given list of Motifs and adds to output.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.info['SlimCheck']: return
            self.obj['SlimCheck'] = rje_slimlist.SLiMList(self.log,self.cmd_list+['motifs=%s' % self.info['SlimCheck']])
            self.obj['SlimCheck'].loadMotifs()
            if not self.obj['SlimCheck'].motifs(): return
            if self.info['ResFile'].lower() in ['','none']:
                self.log.errorLog('Cannot check SLiMs without resfile=FILE',printerror=False)
                return

            ### ~ [2] Make into SLiMFinder SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #X#aa = string.split('A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,^,$',',')
            checklist = []
            for Motif in self.obj['SlimCheck'].slims()[0:]:
                if Motif.slim(): checklist.append(Motif.slim())
                else: self.obj['SlimCheck'].removeMotif(Motif)
            if not self.obj['SlimCheck'].slims(): return

            ### ~ [3] Calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            X = self.obj['SlimCheck'].motifNum() * self.UPNum()
            x = 0.0
            self.list['SlimCheckExtra'] = []   # Extra SLiMs
            for Motif in self.obj['SlimCheck'].motifs()[0:]:
                self.log.printLog('\r#CHECK','SLiM Check (%d motifs) %.2f%%' % (self.obj['SlimCheck'].motifNum(),x/X),newline=False,log=False)
                slim = Motif.slim()
                if slim not in self.dict['Slim']:
                    try:
                        ux = 0
                        sx = 0
                        ox = 0
                        for upc in self.list['UP']:
                            uph = False
                            if self.QUP(upc): continue  # Ignore Query UPC
                            for Seq in upc:
                                hits = Motif.searchSequence(sequence=Seq.info['MaskSeq'])
                                if hits:
                                    uph = True
                                    sx += 1
                                    ox += len(hits)
                            if uph: ux += 1
                        self.dict['Slim'][slim] = {'Occ':[1] * ox,'UP':[1] * ux,'Support':sx}
                        self.slimProb(slim)
                        e = self.dict['Extremf.'][slimPos(slim)] * self.dict['Slim'][slim]['Prob']
                        self.dict['Slim'][slim]['Sig'] = rje.poisson(1,e,callobj=self)
                        self.list['SlimCheckExtra'].append(slim)
                    except:
                        self.log.errorLog('SLiM Error (%s)' % slim,quitchoice=False)
                        self.obj['SlimCheck'].removeMotif(Motif)
                        continue
                else:
                    self.dict['Slim'][slim]['Support'] = self.slimSeqNum(slim)
                    self.dict['Slim'][slim]['UP'] = [1] * self.slimUP(slim)
                    if slim not in self.list['SigSlim']: self.list['SlimCheckExtra'].append(slim)
            self.log.printLog('\r#CHECK','SLiM Check (%d motifs) complete: %d extra SLiMs.' % (self.obj['SlimCheck'].motifNum(),len(self.list['SlimCheckExtra'])))
                    
            ### ~ [4] Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] Setup Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.info['ResFile']))
            reshead = self.resHead()
            totalaa = 0     # Total number of AA in dataset
            for seq in self.seqs():
                if self.opt['Masked']: seq.info['Sequence'] = seq.info['MaskSeq'][0:]
                totalaa += seq.nonX()
            masking = self.maskText()   # Summary of Masking Options
            ## ~ [4b] Output for each checks SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add SliMCalc, OccFilter and StatFilter #!#
            for Motif in self.obj['SlimCheck'].motifs():
                slim = Motif.slim()
                pattern = patternFromCode(slim)
                datadict = {'Dataset':self.dataset(),'RunID':self.info['RunID'],
                            'Masking':masking,'Build':self.myPickle(),'RunTime':self.info['RunTime'],#x#'Date':self.info['Date'],
                            'SeqNum':self.seqNum(),'UPNum':self.UPNum(),'AANum':totalaa,
                            'MotNum':len(self.dict['Slim']),'Rank':'*',
                            'Pattern':pattern,'Occ':len(self.dict['Slim'][slim]['Occ']),'Support':self.slimSeqNum(slim),
                            'IC':self.slimIC(slim),'UP':self.slimUP(slim), 'Norm':self.slimUP(slim),
                            }
                if self.dict['Clouds'] and self.dict['Slim'][slim].has_key('Cloud'):
                    try:
                        datadict['Cloud'] = self.dict['Slim'][slim]['Cloud']
                        datadict['CloudSeq'] = len(self.dict['Clouds'][self.dict['Slim'][slim]['Cloud']]['Seq'])
                        datadict['CloudUP'] = len(self.dict['Clouds'][self.dict['Slim'][slim]['Cloud']]['UPC'])
                    except: pass
                for p in ['ExpUP','Prob','Sig','E']:
                    if self.dict['Slim'][slim].has_key(p): datadict[p] = rje_slim.expectString(self.dict['Slim'][slim][p])
                for h in self.list['Headers']:
                    if Motif and h not in datadict and h in Motif.stat: datadict[h] = Motif.stat[h]
                rje.delimitedFileOutput(self,self.info['ResFile'],reshead,delimit,datadict)

            self.log.printLog('\r#OUT','SLiMCheck Output for %d motifs into %s complete.' % (self.obj['SlimCheck'].motifNum(),self.info['ResFile']))

        except: self.log.errorLog('Major problem with SLiMFinder.slimCheck(%s)' % self.info['SlimCheck'])
#########################################################################################################################
### End of SECTION II: SLiMFinder Class                                                                                 #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def patternFromCode(slim): return rje_slim.patternFromCode(slim)  ### Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG)
#########################################################################################################################
def slimPos(slim): return (string.count(slim,'-') / 2) + 1  ### Returns the number of positions in a slim
#########################################################################################################################
def slimLen(slim): return len(patternFromCode(slim))    ### Returns length of slim
#########################################################################################################################
def slimDif(slim1,slim2): return rje_slimcore.slimDif(slim1,slim2)  ### Returns no. of dif. pos. between slim1 and slim2
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
    try: QSLiMFinder(mainlog,cmd_list).run()
        
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
