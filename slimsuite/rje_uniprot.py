#!/usr/bin/python

# rje_uniprot.py - RJE Module to Handle Uniprot Files
# Copyright (C) 2006 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_uniprot
Description:  RJE Module to Handle Uniprot Files
Version:      3.8
Last Edit:    15/11/10
Copyright (C) 2007 Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains methods for handling UniProt files, primarily in other rje modules but also with some
    standalone functionality. To get the most out of the module with big UniProt files (such as the downloads from EBI),
    first index the UniProt data using the rje_dbase module.

    This module can be used to extract a list of UniProt entries from a larger database and/or to produce summary tables
    from UniProt flat files.

    In addition to method associated with the classes of this module, there are a number of methods that are called from
    the rje_dbase module (primarily) to download and process the UniProt sequence database.

Input Options:
    unipath=PATH    : Path to UniProt Datafile (will look here for DB Index file made with rje_dbase)
    dbindex=FILE    : Database index file [uniprot.index]
    uniprot=FILE    : Name of UniProt file [None]
    extract=LIST    : Extract IDs/AccNums in list. LIST can be FILE or list of IDs/AccNums X,Y,.. []
    acclist=LIST    : As Extract.
    specdat=LIST    : Make a UniProt DAT file of the listed species from the index (over-rules extract=LIST) []
    splicevar=T/F   : Whether to search for AccNum allowing for splice variants (AccNum-X) [True]
    tmconvert=T/F   : Whether to convert TOPO_DOM features, using first description word as Type [False]

Output Options:
    makeindex=T/F   : Generate UniProt index files [False]
    makespec=T/F    : Generate species table [False]
    makefas=T/F     : Generate fasta files [False]
    datout=FILE     : Name of new (reduced) UniProt DAT file of extracted sequences [None]
    tabout=FILE     : Table of extracted UniProt details [None]
    linkout=FILE    : Table of extracted Database links [None]
    longlink=T/F    : Whether link table is to be "long" (acc,db,dbacc) or "wide" (acc, dblinks) [True]
    ftout=FILE      : Outputs table of features into FILE [None]
    domtable=T/F    : Makes a table of domains from uniprot file [False]
    cc2ft=T/F       : Extra whole-length features added for TISSUE and LOCATION (not in datout) [False]

UniProt Conversion Options:
    ucft=X          : Feature to add for UpperCase portions of sequence []
    lcft=X          : Feature to add for LowerCase portions of sequence []
    maskft=LIST     : List of Features to mask out []
    invmask=T/F     : Whether to invert the masking and only retain maskft features [False]
    caseft=LIST     : List of Features to make upper case with rest of sequence lower case []

General Options:
    append=T/F      : Append to results files rather than overwrite [False]
    memsaver=T/F    : Memsaver option to save memory usage - does not retain entries in UniProt object [False]
    cleardata=T/F   : Whether to clear unprocessed Entry data (True) or (False) retain in Entry & Sequence objects [True]

Uses general modules: glob, os, re, string, sys, time
Uses RJE modules: rje, rje_sequence
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, re, string, sys, time
#########################################################################################################################
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_sequence, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial working version for interaction_motifs.py
    # 1.1 - Minor tidying and modification
    # 2.0 - Moved functions to rje_dbase. Added option to extract using index files.
    # 2.1 - Added possibility to extract splice variants
    # 2.2 - Added output of feature table for the entries in memory (not compatible with memsaver mode)
    # 2.3 - Added ID to tabout and also added accShortName() method to extract dictionary of {acc:ID__PrimaryAcc}
    # 2.4 - Add method for converting Sequence object and dictionary into UniProt objects... and saving
    # 2.5 - Added cc2ft Extra whole-length features added for TISSUE and LOCATION [False] and ftout=FILE
    # 2.6 - Added features based on case of sequence. (Uses seq.dict['Case'])
    # 2.7 - Added masking of features - Entry.maskFT(type='EM',inverse=False)
    # 2.8 - Added making of Taxa-specific databases using a list of UniProt Species codes
    # 2.9 - Added extraction of EnsEMBL, HGNC, UniProt and EntrezGene from IPI DAT file.
    # 3.0 - Added some module-level methods for use with rje_dbase.
    # 3.1 - Added extra linking of databases from UniProt entries
    # 3.2 - Added feature masking and TM conversion.
    # 3.3 - Added DBase processing options.
    # 3.4 - Made modifications to allow extended EMBL functionality as part of rje_embl.
    # 3.5 - Added SplitOut to go with rje_embl V0.1
    # 3.6 - Added longlink=T/F  : Whether link table is to be "long" (acc,db,dbacc) or "wide" (acc, dblinks) [True]
    # 3.7 - Added cleardata=T/F : Whether to clear unprocessed Entry data or retain in Entry & Sequence objects [True]
    # 3.8 - Added extraction of NCBI Taxa ID.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Lots of functionality to add! Look also to BioPython.
    # [Y] : Modify the searching for entry in acclist to cope with partial matches (exclude them)
    # [ ] : Modify DomTable to work with Memsaver
    # [ ] : Add specific database detail extraction, first for IPI DAT files and later for UniProt
    # [ ] : Add a database mapping method for extracting DB cross-refs.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_UNIPROT', '3.7', 'January 2010', '2007')
    description = 'RJE Uniprot Parsing/Extraction Module'
    author = 'Dr Richard J. Edwards.'
    comments = []
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: MODULE CONSTANTS                                                                                        #
#########################################################################################################################
### UniProt Parsing dictionary: Add crucial information to parse out here. Used by UniProtEntry.process()   ###
uniparse = {
    'ID' : string.join(['^(\S+)','(\S.+);','(\d+)\s+(AA|BP)'], '\s+'),    # ID, Type, Length
    'AC' : '(\S+);',     # Primary Acc
    'DE' : '\s*(\S.+)',  # Description
    'GN' : 'Name=(\S+);',   # Gene Name
    'SY' : 'Synonyms=(\S.+);',   # Gene Synonyms
    'OS' : '^(\S.+)\s*$',   # Species
    'OX' : '^OX\s+NCBI_TaxID=(\d+)',    # NCBI Taxa ID
    'RX' : 'PubMed=(\d+);', # PubMed ID
    'RC' : 'TISSUE=(.+);',  # Tissue(s)
    'DR' : '^(\S+);\s+(\S.+)$',  # Database links (Dbase,Details)
    'CC' : '^-!-\s+(\S.+):\s+(\S.+)$',  # Comments (Type, Details)
    'FT' : string.join(['(\S+)','<*(\d+)','>*(\d+)\.*','(\S.+)\s*$'], '\s+')   # Feature type, start, stop, desc
    }
#########################################################################################################################
useful_data = ['ID','AC','DE','GN','OS','OC','OX','RX','CC','DR','RC','KW','FT']     # Data to retain following parsing # ?? #
#!# NB. This list is not currently used! #!#
#########################################################################################################################
featurelist = ['LIPID','TRANSMEM','MOD_RES','DOMAIN']   #!# Features for function table. Add more! #!#
#########################################################################################################################
### END OF SECTION II                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: UniProt Class                                                                                          # 
#########################################################################################################################
class UniProt(rje.RJE_Object):     
    '''
    UniProt Download Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of UniProt File 
    - UniPath = Path to UniProt Datafile (will look here for DB Index file) [UniProt/]
    - DBIndex = Database index file [uniprot.index]
    - DATOut = Name of new (reduced) UniProt DAT file of extracted sequences [None]
    - TabOut = Name of table of extracted sequence details [None]
    - LinkOut = Table of extracted Database links [None]
    - FTOut = Outputs table of features into FILE [None]
    - SplitOut = If path given, will split output into individual files per entry into PATH []
    - UCFT = Feature to add for UpperCase portions of sequence []
    - LCFT = Feature to add for LowerCase portions of sequence []
    
    Opt:boolean
    - ClearData = Whether to clear unprocessed Entry data (True) or (False) retain in Entry & Sequence objects [True]
    - DomTable = Makes a table of domains from uniprot file [False]
    - LongLink = Whether link table is to be "long" (acc,db,dbacc) or "wide" (acc, dblinks) [True]    
    - MakeIndex = Generate UniProt index files [False]
    - MakeSpec = Generate species table [False]
    - MakeFas = Generate fasta files [False]
    - SpliceVar = Whether to search for AccNum allowing for splice variants (AccNum-X) [True]

    Stat:numeric

    List:list
    - Entry = list of UniProt Entries
    - Extract = Extract AccNums/IDs in list. LIST can be FILE or list of AccNums X,Y,.. []
    - SpecDat = Make a UniProt DAT file of the listed species from the index []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
    ### Attributes
    def entryNum(self): return len(self.list['Entry'])
#####################4####################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','UniPath','DBIndex','DATOut','TabOut','LinkOut','FTOut','UCFT','LCFT','SplitOut']
        self.optlist = ['DomTable','LongLink','MakeIndex','MakeSpec','MakeFas','SpliceVar','ClearData']
        self.statlist = []
        self.listlist = ['Extract','Entry','SpecDat']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'UniPath':rje.makePath('UniProt/'),'DBIndex':'uniprot.index','UCFT':'','LCFT':''})
        self.setOpt({'SpliceVar':True,'LongLink':True,'ClearData':True})
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
                ### Class Options ###
                self._cmdReadList(cmd,'path',['UniPath'])
                self._cmdReadList(cmd,'file',['DBIndex','DATOut','TabOut','LinkOut','FTOut'])
                self._cmdReadList(cmd,'info',['UCFT','LCFT'])
                self._cmdReadList(cmd,'opt',['DomTable','LongLink','MakeIndex','SpliceVar','MakeSpec','MakeFas','ClearData'])
                self._cmdReadList(cmd,'list',['Extract','SpecDat'])
                self._cmdRead(cmd,type='file',att='Name',arg='uniprot')
                self._cmdRead(cmd,type='list',att='Extract',arg='acclist')
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def run(self):  ### Main Run Method if called direct from commandline
        '''Main Run Method if called direct from commandline. Returns True if no Errors, else False.'''
        try:### ~ [1] ~ Make Index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['MakeIndex'] or self.opt['MakeSpec'] or self.opt['MakeFas']:
                processUniProt(self,makeindex=self.opt['MakeIndex'],makespec=self.opt['MakeSpec'],makefas=self.opt['MakeFas'])
                return
                
            ### Get Extract List from species ###
            if self.list['SpecDat']:
                self.extractSpecies()
                
            ### Check for input needs ###
            if self.info['Name'].lower() in ['','none'] and not self.list['Extract']:   # No AccNum!
                self.log.errorLog('No input file or acclist given. Use "help" option for parameters.',printerror=False)
                return False
            
            ### Setup Output Files ###
            for file in ['DATOut','TabOut','LinkOut','FTOut']:
                if self.info[file].lower() == 'none':
                    self.info[file] = ''
                if self.info[file]:
                    rje.backup(self,file)

            ### Extracted details & MemSaver ###
            if (self.info['TabOut'] or self.info['LinkOut'] or self.info['FTOut']) and self.opt['MemSaver']:
                memtext = 'TabOut, LinkOut and FTOut will not function with MemSaver mode.'
                if self.stat['Interactive'] >= 0 and rje.yesNo('%s Switch Memsaver off?' % memtext):
                    self.opt['MemSaver'] = False
                    self.cmd_list += ['memsaver=F']
                else:
                    self.log.printLog('#MEM',memtext)
            
            ### Read UniProt File ###
            self.readUniProt()
            
            ### Special Features ###
            self.tableOutput() #!# Special Temp #!#
            if self.opt['DomTable']:
                if self.opt['MemSaver']:
                    self.log.errorLog('Old domTable() method in operation. Does not work with memsaver=T.',printerror=False)
                else:
                    self.domTable()
            if self.info['FTOut'] and not self.opt['MemSaver']:
                self.ftTable(self.info['FTOut'])
            return True
        except:
            self.log.errorLog('Fundamental error during UniProt.run()')
            return False
#########################################################################################################################
    def extractSpecies(self):   ### Uses index file to convert species codes into list of Accession Numbers
        '''Uses index file to convert species codes into list of Accession Numbers.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            indexfile = self.info['UniPath'] + self.info['DBIndex']
            if not os.path.exists(indexfile):
                self.log.errorLog('Index file "%s" missing. Cannot make Taxa DAT file!' % indexfile,printerror=False)
                sys.exit()
            self.list['Extract'] = []
            ### ~ [1] ~ Read IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            line = INDEX.readline()
            self.progLog('\r#ACC','Extracting IDs for %d taxa: %s found.' % (len(self.list['SpecDat']),rje.integerString(len(self.list['Extract']))))
            while line:
                if rje.matchExp('(\S+_(%s));' % string.join(self.list['SpecDat'],'|'),line):
                    self.list['Extract'].append(rje.matchExp('(\S+_(%s));' % string.join(self.list['SpecDat'],'|'),line)[0])
                    pc = 100.0 * INDEX.tell() / end_pos
                    self.progLog('\r#ACC','Extracting IDs for %d taxa %.1f%%: %s found.' % (len(self.list['SpecDat']),pc,rje.integerString(len(self.list['Extract']))))
                line = INDEX.readline()
            INDEX.close()                    
            self.printLog('\r#ACC','Extracting IDs for %d taxa complete: %s found.' % (len(self.list['SpecDat']),rje.integerString(len(self.list['Extract']))))
        except: self.errorLog('Problem with %s.extractSpecies()' % self,quitchoice=False); raise
#########################################################################################################################
    ### <2> ### Reading UniProt Entry
#########################################################################################################################
    def readUniProt(self,filename=None,clear=True,acclist=[],logft=False,use_index=True,cleardata=None,reformat=False):    ### Reads UniProt download into UniProtEntry objects
        '''
        Reads UniProt download into UniProtEntry objects.
        >> filename:str = UniProt download filename [None]
        >> clear:boolean = Whether to clear self.list['Entry'] before reading [True]
        >> acclist:list of str objects = UniProt accnum or id list to read
        >> logft:boolean [False] = whether to write number of features to log file
        >> use_index:boolean [True] = whether to use index file if present
        >> cleardata:boolean [True] = whether to clear processed data to save memory 
        >> reformat:boolean = whether to save to DatOut using cut-down method
        << True if success, False if fail
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.opt['ClearData']
            ## ~ [0a] ~ Index and DAT Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            indexfile = ''
            if self.info['DBIndex'] != 'None' and use_index: indexfile = self.info['UniPath'] + self.info['DBIndex']
            ## ~ [0b] ~ Alternative File Name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not filename: filename = self.info['Name']
            if rje.checkForFile(filename) and indexfile and not os.path.exists(indexfile):
                self.log.printLog('#DB','Index file %s not found but DAT file given: ignoring index.' % indexfile)
                indexfile = ''
            if not rje.checkForFile(filename) and (not indexfile or not os.path.exists(indexfile)):
                self.log.printLog('#ERR','UniProt file "%s" and index file "%s" not found!' % (filename,indexfile))
                return False
            ## ~ [0c] ~ Setup Extract list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not acclist: acclist = self.list['Extract'][0:]  
            acclist.sort()
            if not acclist and rje.checkForFile(filename): indexfile = ''   # Process whole file
            ## ~ [0d] ~ DatFiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            datfiles = []
            if indexfile:
                datfiles = glob.glob('%s*.dat' % self.info['UniPath']) + glob.glob('%s*.DAT' % self.info['UniPath'])
                if not datfiles:
                    self.printLog('#ERR','No *.dat files in "%s"! (Cannot use index)' % (self.info['UniPath']))
                    indexfile = ''
            if clear: self.list['Entry'] = []

            ### ~ [1] ~ Process whole file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not indexfile: return self._processWholeFile(filename,logft=logft,reformat=reformat)

            ### ~ [2] ~ Process from Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Index Dictionaries Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dat_keys = {}   # Index key:datfile
            dat_dict = {}   # Dictionary of (key:{list of positions:list of accs/ids})
            acc_dict = {}   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            start_pos = 0
            ## ~ [2b] ~ Read and locate DataFiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            line = INDEX.readline()
            while line and rje.matchExp('^#(\d+)=(\S+)',line):
                self.verbose(1,3,line,0)
                start_pos = INDEX.tell()
                data = rje.matchExp('^#(\d+)=(\S+)',line)
                dat_keys[data[0]] = self.info['UniPath'] + data[1]
                dat_dict[data[0]] = {}
                if not os.path.exists(dat_keys[data[0]]):
                    self.errorLog('%s missing. May not extract all sequences.' % dat_keys[data[0]],printerror=False)
                    dat_keys[data[0]] = None
                line = INDEX.readline()
            ## ~ [2c] Build index dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            re_index = '^(\S+);(\S+):(\S+)'
            splicevar = []
            missing = []
            true_start = start_pos
            for acc in acclist:
                ## Search index ##
                ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   
                if ipos < 0:    # Could not find target!
                    if self.opt['SpliceVar'] and rje.matchExp('(\S+)-(\d+)$',acc): splicevar.append(rje.matchExp('(\S+)-(\d+)$',acc)[0])
                    else: missing.append(acc)
                    continue
                ## Update ##
                (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                    dat_dict[key][pos].append(acc)
                else: dat_dict[key][pos] = [acc]
                if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                    acc_dict[acc].append({key:pos})  
                else: acc_dict[acc] = [{key:pos}]   
                self.progLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
            self.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
            ## ~ [2d] ~ Splice Variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['SpliceVar'] and splicevar:
                ## Setup search ##
                splicevar.sort()    
                if self.list['Extract']: self.list['Extract'] += splicevar
                self.log.printLog('#VAR','Looking for %s potential splice variants.' % (rje.integerString(len(splicevar))))
                start_pos = true_start
                for acc in splicevar[0:]:
                    ## Search index ##
                    ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                    if ipos < 0:    # Could not find target!
                        missing.append(acc)
                        continue
                    ## Update ##
                    (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                    (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                    if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                        dat_dict[key][pos].append(acc)
                    else: dat_dict[key][pos] = [acc]
                    if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                        acc_dict[acc].append({key:pos})  
                    else: acc_dict[acc] = [{key:pos}]   
                    self.progLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
                self.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acclist)),rje.integerString(len(missing))))
            ## ~ [2e] ~ Missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for acc in missing: self.printLog('#ACC','AccNum/ID "%s" missing from %s' % (acc,indexfile))
            
            ### ~ [3] ~ Extract From DAT Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            extract_dict = {}   # {{key:pos}:'ID (AccNum)'} = For matching input acclist to extracted entries
            for key in rje.sortKeys(dat_keys):
                unifile = dat_keys[key]
                if not rje.checkForFile(unifile): self.errorLog('%s entries skipped due to missing file' % unifile,printerror=False); continue     
                UNIFILE = open(unifile,'r')
                ex = 0
                for pos in dat_dict[key]:
                    ipos = string.atol(pos)
                    # Get to correct position
                    rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[1]
                    # Read
                    if self._readSingleEntry(UNIFILE,logft=logft,cleardata=cleardata,reformat=reformat):
                        # Check for replaced AccNum/ID
                        wanted_acc = dat_dict[key][pos]
                        new_id = self.list['Entry'][-1].obj['Sequence'].info['ID']
                        new_acc = self.list['Entry'][-1].obj['Sequence'].info['AccNum']
                        for acc in wanted_acc:
                            if acc not in [new_id,new_acc]:
                                if self.stat['Verbose'] > 0: self.printLog('\n#ACC','Secondary AccNum %s mapped to %s (%s).' % (acc,new_id,new_acc))
                                else: self.printLog('#ACC','Secondary AccNum %s mapped to %s (%s).' % (acc,new_id,new_acc),screen=False)
                        # Update
                        ex += 1
                        if self.opt['MemSaver']: self.list['Entry'] = []     # Delete for memsaver
                    else:
                        bummer = rje.chomp(rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[0])
                        self.errorLog('%s rejected by _readSingleEntry() but explicitly selected for extraction!' % bummer,printerror=False)
                    self.progLog('\r#DAT','%s entries extracted from %s.' % (rje.integerString(ex),unifile))
                UNIFILE.close()
                self.printLog('\r#DAT','%s entries extracted from %s.' % (rje.integerString(ex),unifile))
            return True
        except: self.errorLog('UniProt.readUniProt() Failed. Check format.'); return False
#########################################################################################################################
    def _processWholeFile(self,filename,logft=True,cleardata=None,reformat=False,log=True):    ### Processes whole file into entries
        '''
        Processes whole file into entries.
        >> filename:str = UniProt filename
        >> logft:boolean = whether to write number of features to log file
        >> reformat:boolean = whether to save to DatOut using cut-down method
        << returns True/False dependent on success
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.opt['ClearData']
            rx = 0      # Number of entries read
            ex = 0      # Number of entries extracted
            UNIPROT = open(filename, 'r')
            ### ~ [1] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Process'
            logtext = 'Extracting entries from %s' % string.split(filename,os.sep)[-1]
            while self._readSingleEntry(UNIPROT,logft=logft,cleardata=cleardata,reformat=reformat):
                rx += 1
                if self.opt['MemSaver'] and self.entryNum() > 0:    # Kept sequence
                    ex += 1
                    self.list['Entry'] = []     # Delete for memsaver
                else: ex = self.entryNum()
                if log: self.progLog('\r#DAT','%s: %s read, %s extracted.' % (logtext,rje.integerString(rx),rje.integerString(ex)))
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log: self.printLog('\r#DAT','%s: %s read, %s extracted.' % (logtext,rje.integerString(rx),rje.integerString(ex)))    
            return True     
        except: self.errorLog('Cataclysmic error during %s._processWholeFile()!' % self); return False
#########################################################################################################################
    def _newEntryObject(self): return UniProtEntry(log=self.log,cmd_list=self.cmd_list)
#########################################################################################################################
    def _readSingleEntry(self,UNIPROT,logft=True,cleardata=None,reformat=False):    ### Reads a single entry from current point in file and processes 
        '''
        Processes whole file into entries.
        >> UNIPROT:FileHandle = Open UniProt file for reading *at start of entry*
        >> logft:boolean = whether to write number of features to log file
        >> cleardata:boolean = whether to clear data to save memory after processing
        >> reformat:boolean = whether to save to DatOut using cut-down method
        << returns True/False dependent on success
        '''
        try:### ~ [0] ~ Check for end of file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cleardata == None: cleardata = self.opt['ClearData']
            line = UNIPROT.readline()
            if not line: return False
            ### ~ [1] ~ Process Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _entry = None       # Current entry
            _reading = False    # Whether currently reading an entry
            while line:
                ## ~ [1a] ~ Full Text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if _entry: _entry.info['FullText'] = '%s%s' % (_entry.info['FullText'],line)
                line = rje.chomp(line)
                type = line[0:2]
                if type in ['','XX']: line = UNIPROT.readline(); continue
                rest = line[5:]
                ## ~ [1b] ~ New Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                    
                if type == 'ID':
                    _reading = True
                    _entry = self._newEntryObject() 
                    _entry.info['FullText'] = '%s\n' % line
                elif not _entry:
                    self.errorLog('Expected ID entry, got "%s"' % line,printerror=False)
                    raise ValueError
                ## ~ [1c] ~ End of Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                _stage = 'End of Entry'
                if type == '//':
                    if self._add_entry(_entry,acclist=self.list['Extract']):     # Entry is one of the desired entries
                        self.list['Entry'].append(_entry)   # Add entry to list
                        _entry.process(logft=logft,cleardata=cleardata and not reformat)         # Extracts details from uniprot
                        if self.info['DATOut'] and self.info['DATOut'].lower() != 'none':
                            if self.info['SplitOut'].lower() not in ['','none']:
                                datapp = False
                                datout = rje.makePath(self.info['SplitOut']) + self.info['DATOut']
                                rje.mkDir(self,datout)
                                datout = string.join([datout,_entry.info['Type'],_entry.info['ID'],'dat'],'.')
                            else: datout = self.info['DATOut']; datapp = True                            
                            if reformat: self.saveUniProt(datout,[_entry],append=datapp)
                            else: open(datout,{True:'a',False:'w'}[datapp]).write(_entry.info['FullText'])
                        #!# >>>>>> This is a fudge but it's OK for now <<<<<<<<<<< #!#
                        #X#self.tableOutput(_entry)
                        for k in ['DBLinks']:
                            if not self.list.has_key(k): self.list[k] = []
                        for k in ['DBLinks']:
                            for i in _entry.dict[k].keys():
                                if i not in self.list[k]: self.list[k].append(i)
                        #!# ^^^^^^^^^^ This is a fudge but it's OK for now ^^^^^^^ #!#
                        _entry.info['FullText'] = ''
                        return True
                    else: self.deBug('Bad Entry'); self.deBug(_entry); self.deBug(self.list['Extract']); return False
                ## ~ [1d] ~ Entry Details read into a dictionary within the entry ~~~~~~~~~~~~~~~~~ ##
                if not rest: self.deBug(line)
                if _entry.dict['Data'].has_key(type):   # Append list
                    if rest[:1] != ' ': _entry.dict['Data'][type].append(rest)   # New entry
                    else:
                        while rest[:1] == ' ': rest = rest[1:]
                        _entry.dict['Data'][type][-1] = '%s %s' % (_entry.dict['Data'][type][-1], rest)
                elif type == '  ':
                    if _entry.dict['Data'].has_key('SEQ'):
                        _entry.dict['Data']['SEQ'][0] = '%s %s' % (_entry.dict['Data']['SEQ'][0],rest)
                    elif _entry.dict['Data'].has_key('SQ'): _entry.dict['Data']['SEQ'] = [rest]
                elif type not in ['XX']: _entry.dict['Data'][type] = [rest]
                ## ~ [1e] ~ Next Line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                line = UNIPROT.readline()
            ### ~ [2] ~ Reached EOF before // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if _reading: self.errorLog('Started UniProt Entry but EOF reached before "//". Truncated input file?',printerror=False)
            return False                
        except: self.errorLog('Cataclysmic error during %s._readSingleEntry()!' % self,quitchoice=False); raise
#########################################################################################################################
    def _add_entry(self,_entry,acclist):    ### Returns True if _entry in acclist or false if not
        '''
        Returns True if _entry in acclist or false if not.
        >> _entry:uniProtEntry object (unprocessed)
        >> acclist:list of accession numbers and/or IDs
        '''
        try:### ~ Check for entry in acclist ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not acclist: return True
            for acc in acclist:
                if _entry.dict['Data'].has_key('AC') and string.join([' '] + _entry.dict['Data']['AC'],' ').find(' %s;' % acc) > 0: return True
                if _entry.dict['Data'].has_key('ID') and string.join(_entry.dict['Data']['ID'],' ').find('%s ' % acc) == 0: return True
        except: self.errorLog('Cataclysmic error during _add_entry')
        return False
#########################################################################################################################
    def addFromSeq(self,seq=None,sequence='',name='',data={},ft=[]):    ### Converts into UniProtEntry object 
        '''
        Converts into UniProtEntry object and adds to self.
        >> seq:rje_sequence.Sequence object [None]
        >> sequence:str = alternative sequence data (will be converted to Sequence object!) ['']
        >> name:str = alternative sequence name (will be converted to Sequence object!) ['']
        >> data:dict = dictionary of UniProt data with {keys ID/AC/OS etc: [list of lines]} [{}]
        >> ft:list = list of ftdic dictionaries of features {'Type/Desc':str,'Start/End':int} [[]]
        << returns entry if successful or None if fails
        '''
        ### ~ [1] ~ Add case features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        addft = ft[0:]
        if self.info['UCFT'] and seq:
            for uc in seq.dict['Case']['Upper']:
                addft.append({'Type':self.info['UCFT'].upper(),'Desc':self.info['UCFT'],'Start':uc[0]+1,'End':uc[1]+1})
        if self.info['LCFT'] and seq:
            for lc in seq.dict['Case']['Lower']:
                addft.append({'Type':self.info['LCFT'].upper(),'Desc':self.info['LCFT'],'Start':lc[0]+1,'End':lc[1]+1})
        ### ~ [2] ~ Make Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###        
        newentry = self._newEntryObject()
        entry = newentry.uniProtFromSeq(seq,sequence,name,data,addft)
        if entry: self.list['Entry'].append(entry)
        return entry
#########################################################################################################################
    ### <3> ### UniProt Info Output
#########################################################################################################################
    def tableOutput(self):   ### Tabulated output of UniProt information
        '''Tabulated output of UniProt information. Divided into TabOut (UniProt summary) and LinkOut (database links)'''
        try:
            ### Database Links Table ###
            if self.info['LinkOut'].lower() not in ['','none']:
                self.linkOutput()
            if self.info['TabOut'].lower() in ['','none']:
                return

            ### Setup Groups and Columns ###
            headers = {'#1# Basic Details #':['#','AccNum','ID','Len','Description','Gene','Species'],
                       '#2# Function & Activity #':['Function','GO_MF','GO_BP','Activity','Interactions','Phenotype',
                                                    'Similarity'],
                       '#3# Expression, Location & Structure #':['Tissue','Cell_Loc','PDB','InterPro','Pfam','PROSITE',
                                                                 'Isoforms','Ensembl'],
                       '#4# References & Links #':['GeneCards','PubMed','Keywords','Comments']
                       }

            ### Open File and write header ###
            delimit = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['TabOut']))
            TABOUT = open(self.info['TabOut'],'a')   # Already deleted if append=F
            TABOUT.write('# Generated by %s: %s\n' % (self.log.info['Name'],time.asctime(time.localtime(time.time()))))
            if self.info['Name'] != 'None':
                TABOUT.write('# Source: %s\n' % os.path.abspath(self.info['Name']))
            else:
                TABOUT.write('#Source: %s\n' % os.path.abspath(self.info['UniPath']+'*.dat'))
            TABOUT.write('# Seqnum: %s\n\n' % rje.integerString(self.entryNum()))
            ## Headers ##
            head1 = []
            head2 = []
            for h in rje.sortKeys(headers):
                head1.append(h)
                head1 += [''] * (len(headers[h])-1)
                head2 += headers[h]
            rje.writeDelimit(TABOUT,head1,delimit)
            rje.writeDelimit(TABOUT,head2,delimit)

            ### Write Data for Entries ###
            ex = 0
            for entry in self.list['Entry']:
                seq = entry.obj['Sequence']
                ex += 1
                data = []
                comments = rje.sortKeys(entry.dict['Comments'])
                for h in head2:     # Column headers:
                    #1# Basic Details #
                    if h == '#':
                        data.append(rje.preZero(ex,self.entryNum()))
                    elif h in ['AccNum','ID','Description']:
                        data.append(seq.info[h])
                    elif h == 'Len':
                        data.append('%d' % seq.aaLen())
                    elif h == 'Gene':
                        data.append(string.join([seq.info['Gene']] + entry.list['Synonyms'],'; '))
                    elif h == 'Species':
                        data.append('%s [%s]' % (seq.info['Species'],seq.info['SpecCode']))
                    #2# Function & Activity #
                    elif h == 'Function':
                        text = ''
                        for cc in ['FUNCTION','PATHWAY','DOMAIN']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'GO_MF':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; F:') < 0:
                                    go.remove(g)
                        data.append(string.join(go,' >> '))
                    elif h == 'GO_BP':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; P:') < 0:
                                    go.remove(g)
                        data.append(string.join(go,' >> '))
                    elif h == 'Activity':   #!# Join to Function #!#
                        text = ''
                        for cc in ['CATALYTIC ACTIVITY','COFACTOR']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Interactions':
                        text = ''
                        for cc in ['INTERACTION','ENZYME REGULATION','SUBUNIT','PTM']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Phenotype':
                        text = ''
                        for cc in ['DISEASE','POLYMORPHISM']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Similarity':
                        text = ''
                        for cc in ['SIMILARITY']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    #3# Expression, Location & Structure #
                    elif h == 'Tissue':
                        text = ''
                        if entry.list['Tissues']:
                            text = 'TISSUES: %s' % string.join(entry.list['Tissues']+[''],'; ')
                        for cc in ['TISSUE SPECIFICITY','DEVELOPMENTAL STAGE','INDUCTION']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    elif h == 'Cell_Loc':
                        go = []
                        if entry.dict['DBLinks'].has_key('GO'):
                            go = entry.dict['DBLinks']['GO'][0:]
                            for g in go[0:]:
                                if g.find('; C:') < 0:
                                    go.remove(g)
                        cc = 'SUBCELLULAR LOCATION'
                        if entry.dict['Comments'].has_key(cc):
                            comments.remove(cc)
                            go = entry.dict['Comments'][cc] + go
                        data.append(string.join(go,' >> '))
                    elif h in ['PDB','InterPro','Pfam','PROSITE','Ensembl']:
                        if entry.dict['DBLinks'].has_key(h):
                            data.append(string.join(entry.dict['DBLinks'][h],' >> '))
                        else:
                            data.append('')
                    elif h == 'Isoforms':
                        text = ''
                        for cc in ['ALTERNATIVE PRODUCTS']:
                            if entry.dict['Comments'].has_key(cc):
                                comments.remove(cc)
                                if text:
                                    text += ' '
                                text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                                if text[-1] != '.':
                                    text += '.'
                        data.append(text)
                    #4# References & Links
                    elif h == 'GeneCards':
                        if entry.dict['DBLinks'].has_key('HGNC'):
                            data.append(string.join(entry.dict['DBLinks']['HGNC'],' >> '))
                        else:
                            data.append('')
                    elif h in ['PubMed']:
                        if entry.list[h]:
                            data.append('http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=%s' % string.join(entry.list[h],','))
                        else:
                            data.append('')
                    elif h in ['Keywords']:
                        if entry.list[h]:
                            data.append(string.join(entry.list[h],','))
                        else:
                            data.append('')
                    elif h == 'Comments':
                        text = ''
                        for cc in comments:
                            if text:
                                text += ' '
                            text += '%s: %s' % (cc,string.join(entry.dict['Comments'][cc],' >> '))
                            if text[-1] != '.':
                                text += '.'
                        data.append(text)

                rje.writeDelimit(TABOUT,data,delimit)
                self.log.printLog('\r#OUT','UniProt summary output: %.1f%%' % (100.0 * ex / self.entryNum()),log=False,newline=False)
            self.log.printLog('\r#OUT','UniProt summary output for %s entries.' % rje.integerString(self.entryNum()))
            TABOUT.close()


        except:
            self.log.errorLog('Major error during UniProt.tableOutput()!',quitchoice=True)            
#########################################################################################################################
    def linkOutput(self):   ### Delimited output of UniProt database links
        '''Delimited output of UniProt database links.'''
        try:
            ### Setup ##
            if self.opt['LongLink']: self.linkOutputLong()  # self.obj['Sequence'].list['Secondary ID']
            delimit = rje.getDelimit(self.cmd_list,default=rje.delimitFromExt(filename=self.info['LinkOut']))
            dblist = rje.dictValues(self.list,'DBLinks')
            dblist.sort()   #!# Group Databases by type later #!#
            
            ### Open File and write header ###
            LINKFILE = open(self.info['LinkOut'],'a')   # Already deleted if append=F
            LINKFILE.write('# Generated by %s: %s\n' % (self.log.info['Name'],time.asctime(time.localtime(time.time()))))
            if self.info['Name'] != 'None':
                LINKFILE.write('# Source: %s\n' % os.path.abspath(self.info['Name']))
            else:
                LINKFILE.write('# Source: %s\n' % os.path.abspath(self.info['UniPath']+'*.dat'))
            LINKFILE.write('# Seqnum: %s\n\n' % rje.integerString(self.entryNum()))
            rje.writeDelimit(LINKFILE,['#','AccNum']+dblist,delimit)
            
            ### Write Data for Entries ###
            ex = 0
            for entry in self.list['Entry']:
                ex += 1
                data = [rje.preZero(ex,self.entryNum()),entry.obj['Sequence'].info['AccNum']]
                for db in dblist:
                    if entry.dict['DBLinks'].has_key(db):
                        data.append(string.join(entry.dict['DBLinks'][db],' >> '))
                    else:
                        data.append('')
                rje.writeDelimit(LINKFILE,data,delimit)
                self.log.printLog('\r#OUT','Database Links output: %.1f%%' % (100.0 * ex / self.entryNum()),log=False,newline=False)
            self.log.printLog('\r#OUT','Database Links output for %s entries.' % rje.integerString(self.entryNum()))
            LINKFILE.close()
            
        except:
            self.log.errorLog('Major error during UniProt.linkOutput()!',quitchoice=True)            
#########################################################################################################################
    def linkOutputLong(self):   ### Delimited output of UniProt database links
        '''Delimited output of UniProt database links.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            listhead = ['AccNum','DBase','LinkAcc']
            rje.delimitedFileOutput(self,self.info['LinkOut'],listhead,rje_backup=True)
            ### ~ [1] ~ Write Data for Entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0
            for entry in self.list['Entry']:
                ex += 1
                for secid in self.obj['Sequence'].list['Secondary ID']:
                    data = {'AccNum':entry.obj['Sequence'].info['AccNum'],'DBase':'UniProt','LinkAcc':secid}
                    rje.delimitedFileOutput(self,self.info['LinkOut'],listhead,datadict=data)
                for db in rje.sortKeys(entry.dict['DBLinks']):
                    for dbid in entry.dict['DBLinks'][db]:
                        data = {'AccNum':entry.obj['Sequence'].info['AccNum'],'DBase':db,'LinkAcc':dbid}
                        rje.delimitedFileOutput(self,self.info['LinkOut'],listhead,datadict=data)
                self.progLog('\r#OUT','Database Links output: %.1f%%' % (100.0 * ex / self.entryNum()))
            self.printLog('\r#OUT','Database Links output for %s entries.' % rje.integerString(self.entryNum()))
        except: self.errorLog('Major error during UniProt.linkOutput()!',quitchoice=True)            
#########################################################################################################################
    def domTable(self): ### Outputs domain info into a table    #!# OLD and Unchecked #!#
        '''Outputs domain info into a table.'''
        try:
            ### Setup ###
            delimit = rje.getDelimit(self.cmd_list)
            if self.info['Name'].lower() not in ['','none']:
                outfile = '%s.domains.%s' % (rje.baseFile(self.info['Name']),rje.delimitExt(delimit))
            else:
                outfile = '%s.domains.%s' % (rje.baseFile(self.info['TabOut']),rje.delimitExt(delimit))
            if self.opt['Append']:
                DOMFILE = open(outfile,'a')
            else:
                DOMFILE = open(outfile,'w')
                rje.writeDelimit(DOMFILE,['acc_num','domain','dom_start','dom_end'],delimit)

            ### Output ###
            dx = 0
            for entry in self.list['Entry']:
                seq = entry.obj['Sequence']
                acc = seq.info['AccNum']
                for ft in entry.list['Feature']:
                    if ft['Type'].upper() != 'DOMAIN':
                        continue
                    outlist = [acc]
                    for fk in ['Desc','Start','End']:
                        outlist.append('%s' % ft[fk])
                    rje.writeDelimit(DOMFILE,outlist,delimit)
                    dx += 1

            ### End ###
            self.log.printLog('#DOM','%d domains from %d proteins output to %s.' % (dx,len(self.list['Entry']),outfile))
            DOMFILE.close()
                    
        except:
            self.log.errorLog('Cataclysmic error during domTable!')
#########################################################################################################################
    def ftTable(self,outfile):  ### Outputs features into a table
        '''
        Outputs features into a table.
        >> outfile:str = Name of output file
        '''
        try:
            ### Setup ###
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=outfile))
            if self.opt['Append']:
                FT = open(outfile,'a')
            else:
                FT = open(outfile,'w')
                rje.writeDelimit(FT,['acc_num','feature','ft_start','ft_end','description'],delimit)

            ### Output ###
            (fx,ex) = (0,0.0)
            for entry in self.list['Entry']:
                ex += 100.0
                acc = entry.obj['Sequence'].info['AccNum']
                ## Make dictionary of {start:{end:[features]}}
                ft_dict = {}
                for ft in entry.list['Feature']:
                    ft_start = ft['Start']
                    if not ft_dict.has_key(ft_start):
                        ft_dict[ft_start] = {}
                    ft_end = ft['End']
                    if not ft_dict[ft_start].has_key(ft_end):
                        ft_dict[ft_start][ft_end] = []
                    ft_dict[ft_start][ft_end].append(ft)
                ## Sort and output ##
                for ft_start in rje.sortKeys(ft_dict):
                    for ft_end in rje.sortKeys(ft_dict[ft_start]):
                        for ft in ft_dict[ft_start][ft_end]:
                            outlist = [acc]
                            for fk in ['Type','Start','End','Desc']:
                                outlist.append('%s' % ft[fk])
                            rje.writeDelimit(FT,outlist,delimit)
                        fx += 1
                self.log.printLog('\r#FT','Feature output: %.1f%% (%s features)' % (ex/len(self.list['Entry']),rje.integerString(fx)),log=False,newline=False)

            ### End ###
            FT.close()
            self.log.printLog('\r#FT','Feature output complete: %s features, %s entries.' % (rje.integerString(fx),rje.integerString(len(self.list['Entry']))))
        except:
            self.log.errorLog('Program error during rje_uniprot.ftTable()',quitchoice=True)
#########################################################################################################################
    def saveUniProt(self,outfile,entries=[],append=False):    ### Saves self as a DAT file
        '''
        Saves self as a DAT file.
        >> outfile:str = Name of output file
        >> entries:list of entries (self.list['Entry'] if none given)
        >> append:boolean = whether to append file
        '''
        try:### ~ [0] ~ Setup Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not append: rje.backup(self,outfile)
            if not entries: entries = self.list['Entry'][0:]
            ### ~ [1] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OUT = open(outfile,'a')                    
            for entry in entries[0:]:
                if not entry.dict['Data'] and not entry.uniProtFromSeq():
                    entries.remove(entry)
                    self.errorLog('Problem with %s (%s) - cannot output' % (entry,entry.info['Name']),printerror=False)
                    continue
                seq = entry.obj['Sequence']
                ## ~ [1a] ~ Standard info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for key in ['ID','AC','DT','DE','GN','OS']:
                    if entry.dict['Data'].has_key(key):
                        for rest in entry.dict['Data'][key]: OUT.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## ~ [1b] ~ Other data, except Features and sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for key in rje.sortKeys(entry.dict['Data']):
                    if key not in ['ID','AC','DT','DE','GN','OS','FT','SQ','SEQ','//']:
                        for rest in entry.dict['Data'][key]: OUT.write('%s   %s\n' % (key,rje.chomp(rest)))
                ## ~ [1c] ~ Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry.orderFT()
                for ftdict in entry.list['Feature']:
                    (p1,p2) = (ftdict['Start'],ftdict['End'])
                    ftxt = 'FT   %s' % ftdict['Type']
                    while len(ftxt) < 14 or ftxt[-1] != ' ': ftxt += ' '
                    ftxt += '%6s' % ('%d' % p1)
                    while len(ftxt) > 20 and ftxt[-(len('%d' % p1)+2):-len('%d' % p1)] == '  ': ftxt = ftxt[:-(len('%d' % p1)+1)] + ftxt[-len('%d' % p1):]
                    ftxt += '%7s' % ('%d' % p2)
                    while len(ftxt) > 27 and ftxt[-(len('%d' % p2)+2):-len('%d' % p2)] == '  ': ftxt = ftxt[:-(len('%d' % p2)+1)] + ftxt[-len('%d' % p2):]
                    ftxt += ' %s\n' % ftdict['Desc']
                    OUT.write(ftxt)
                ## ~ [1d] ~ Sequence/End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if seq.dna(): OUT.write('SQ   SEQUENCE%s%d BP;  XXX MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % seq.aaLen())),seq.aaLen()))
                else: OUT.write('SQ   SEQUENCE%s%d AA;  %d MW;  000000000000000 RJE06;\n' % (' ' * (7 - len('%d' % seq.aaLen())),seq.aaLen(),rje_sequence.MWt(seq.info['Sequence'])))
                uniseq = seq.info['Sequence'][0:]
                while len(uniseq) > 0:
                    OUT.write('     %s\n' % string.join([uniseq[0:10],uniseq[10:20],uniseq[20:30],uniseq[30:40],uniseq[40:50],uniseq[50:60]],' '))
                    uniseq = uniseq[60:]
                OUT.write('//\n')
            OUT.close()
            if not append or len(entries) > 1: self.printLog('#OUT','UniProt format for %d entries saved to %s' % (len(entries),outfile))            
        except: self.errorLog('Major problem with %s.saveUniProt()' % self)
#########################################################################################################################
    ### <4> ### Additional UniProt Tools                                                                                #
#########################################################################################################################
    def accNameSeq(self,acc_list=[],spec=None,justsequence=True):  ### Method to extract dictionaries of {acc:'ID__PrimaryAcc Desc'} & {acc:seq}
        '''
        Method to extract dictionary of {acc:'ID__PrimaryAcc Desc'} & {acc:seq} using index.
        >> acclist:list of accession numbers. Will use self.list['Extract'] if none given
        >> spec:Limit to species code
        >> justsequence:bool [True] = Whether to just return the sequence data (True) or a sequence object (False)
        << tuple of dictionaries of ({acc:ID__PrimaryAcc},{acc:uniprot sequence})
        '''
        try:
            ### Setup ###
            accshortname = {}
            accseq = {}
            ## Index and DAT Files ##
            indexfile = ''
            if self.info['DBIndex'] != 'None':
                indexfile = self.info['UniPath'] + self.info['DBIndex']
            ## Extract list ##
            if acc_list:
                acc_list.sort() #X#
            else:
                acc_list = self.list['Extract'][0:]  #X# 
                acc_list.sort()
            ## DatFiles ##
            datfiles = glob.glob('%s*.dat' % self.info['UniPath']) + glob.glob('%s*.DAT' % self.info['UniPath'])
            if not datfiles:
                self.log.printLog('#ERR','No *.dat files in "%s"! (Cannot use index)' % (self.info['UniPath']))
                return {}

            ### Process from Index ###
            _stage = 'Index Dictionaries (Setup)'
            dat_keys = {}   # Index key:datfile
            dat_dict = {}   # Dictionary of (key:{list of positions:list of accs/ids})
            acc_dict = {}   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            start_pos = 0
            ## Read DataFiles ##
            line = INDEX.readline()
            while line and rje.matchExp('^#(\d+)=(\S+)',line):
                self.verbose(1,3,line,0)
                start_pos = INDEX.tell()
                data = rje.matchExp('^#(\d+)=(\S+)',line)
                dat_keys[data[0]] = self.info['UniPath'] + data[1]
                dat_dict[data[0]] = {}
                if not os.path.exists(dat_keys[data[0]]):
                    self.log.errorLog('%s missing. May not extract all sequences.' % dat_keys[data[0]],printerror=False)
                    dat_keys[data[0]] = None
                #X#print rje.fileLineFromSeek(INDEX,start_pos,reseek=False,next=False)
                line = INDEX.readline()

            ### Build index dictionary ###
            _stage = 'Build Index Dictionaries'
            re_index = '^(\S+);(\S+):(\S+)'
            splicevar = []
            missing = []
            true_start = start_pos
            for acc in acc_list:
                ## Search index ##
                ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                if ipos < 0:    # Could not find target!
                    if self.opt['SpliceVar'] and rje.matchExp('(\S+)-(\d+)$',acc):
                        splicevar.append(rje.matchExp('(\S+)-(\d+)$',acc)[0])
                    else:
                        missing.append(acc)
                    continue
                ## Update ##
                (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                    dat_dict[key][pos].append(acc)
                else:
                    dat_dict[key][pos] = [acc]
                if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                    acc_dict[acc].append({key:pos})  
                else:                    
                    acc_dict[acc] = [{key:pos}]   
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
            self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Splice Variants ###
            if self.opt['SpliceVar'] and splicevar:
                ## Setup search ##
                splicevar.sort()    #X#                splicevar = rje.sortUnique(splicevar)
                if self.list['Extract']:
                    self.deBug(self.list['Extract'])
                    self.list['Extract'] += splicevar
                    self.deBug(self.list['Extract'])
                self.log.printLog('#VAR','Looking for %s potential splice variants.' % (rje.integerString(len(splicevar))))
                start_pos = true_start
                for acc in splicevar[0:]:
                    ## Search index ##
                    ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                    if ipos < 0:    # Could not find target!
                        missing.append(acc)
                        continue
                    ## Update ##
                    (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                    (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                    if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                        dat_dict[key][pos].append(acc)
                    else:
                        dat_dict[key][pos] = [acc]
                    if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                        acc_dict[acc].append({key:pos})  
                    else:                    
                        acc_dict[acc] = [{key:pos}]   
                    self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Missing acc have no entry in returned dictionary ###
            #X#for acc in missing:            
            #X#    self.log.printLog('#ACC','AccNum/ID "%s" missing from %s' % (acc,indexfile))
            
            ### Extract From UniProt using dictionaries ###
            _stage = 'Extract using Dictionaries'
            extract_dict = {}   # {{key:pos}:'ID (AccNum)'} = For matching input acclist to extracted entries
            for key in rje.sortKeys(dat_keys):
                unifile = dat_keys[key]
                UNIFILE = open(unifile,'r')
                ex = 0
                for pos in dat_dict[key]:
                    ipos = string.atol(pos)
                    # Get to correct position
                    rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[1]
                    # Read
                    if self._readSingleEntry(UNIFILE,logft=False,cleardata=True):
                        # Check for replaced AccNum/ID
                        wanted_acc = dat_dict[key][pos]
                        new_id = self.list['Entry'][-1].obj['Sequence'].info['ID']
                        new_acc = self.list['Entry'][-1].obj['Sequence'].info['AccNum']
                        new_desc = self.list['Entry'][-1].obj['Sequence'].info['Description']
                        if spec and self.list['Entry'][-1].obj['Sequence'].info['SpecCode'] != spec:
                            self.list['Entry'] = self.list['Entry'][:-1]     # Delete 
                            continue                            
                        #X#print new_id, new_acc, wanted_acc
                        for acc in wanted_acc:
                            accshortname[acc] = '%s__%s %s' % (new_id,new_acc,new_desc)
                            if justsequence: accseq[acc] = self.list['Entry'][-1].obj['Sequence'].info['Sequence']
                            else: accseq[acc] = self.list['Entry'][-1].obj['Sequence']
                        # Update
                        ex += 1
                        self.list['Entry'] = self.list['Entry'][:-1]     # Delete 
                    else:
                        bummer = rje.chomp(rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[0])
                        self.log.errorLog('%s rejected by _readSingleEntry() but explicitly selected for extraction!' % bummer,printerror=False)
                    self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile),log=False,newline=False)
                UNIFILE.close()
                self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile))

            return (accshortname,accseq)
            
        except:
            self.log.errorLog('Program error during rje_uniprot.accShortName()',quitchoice=True)
#########################################################################################################################
    def accDict(self,acc_list=[],cleardata=None):      ### Method to extract dictionaries of {acc:UniProtEntry}
        '''
        Method to extract dictionaries of {acc:UniProtEntry}.
        >> acclist:list of accession numbers. Will use self.list['Extract'] if none given
        << dictionary of {acc:UniProtEntry}
        '''
        try:
            ### Setup ###
            if cleardata == None: cleardata = self.opt['ClearData']
            accentry = {}
            ## Index and DAT Files ##
            indexfile = ''
            if self.info['DBIndex'] != 'None': indexfile = self.info['UniPath'] + self.info['DBIndex']
            ## Extract list ##
            if acc_list: acc_list.sort() #X#
            else:
                acc_list = self.list['Extract'][0:]  #X# 
                acc_list.sort()
            ## DatFiles ##
            datfiles = glob.glob('%s*.dat' % self.info['UniPath']) + glob.glob('%s*.DAT' % self.info['UniPath'])
            if not datfiles:
                self.log.printLog('#ERR','No *.dat files in "%s"! (Cannot use index)' % (self.info['UniPath']))
                return {}

            ### Process from Index ###
            _stage = 'Index Dictionaries (Setup)'
            dat_keys = {}   # Index key:datfile
            dat_dict = {}   # Dictionary of (key:{list of positions:list of accs/ids})
            acc_dict = {}   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
            INDEX = open(indexfile,'r')
            INDEX.seek(0,2)
            end_pos = INDEX.tell()
            INDEX.seek(0)
            start_pos = 0
            ## Read DataFiles ##
            line = INDEX.readline()
            while line and rje.matchExp('^#(\d+)=(\S+)',line):
                self.verbose(1,3,line,0)
                start_pos = INDEX.tell()
                data = rje.matchExp('^#(\d+)=(\S+)',line)
                dat_keys[data[0]] = self.info['UniPath'] + data[1]
                dat_dict[data[0]] = {}
                if not os.path.exists(dat_keys[data[0]]):
                    self.log.errorLog('%s missing. May not extract all sequences.' % dat_keys[data[0]],printerror=False)
                    dat_keys[data[0]] = None
                #X#print rje.fileLineFromSeek(INDEX,start_pos,reseek=False,next=False)
                line = INDEX.readline()

            ### Build index dictionary ###
            _stage = 'Build Index Dictionaries'
            re_index = '^(\S+);(\S+):(\S+)'
            splicevar = []
            missing = []
            true_start = start_pos
            for acc in acc_list:
                ## Search index ##
                ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                if ipos < 0:    # Could not find target!
                    if self.opt['SpliceVar'] and rje.matchExp('(\S+)-(\d+)$',acc):
                        splicevar.append(rje.matchExp('(\S+)-(\d+)$',acc)[0])
                    else:
                        missing.append(acc)
                    continue
                ## Update ##
                (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                    dat_dict[key][pos].append(acc)
                else:
                    dat_dict[key][pos] = [acc]
                if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                    acc_dict[acc].append({key:pos})  
                else:                    
                    acc_dict[acc] = [{key:pos}]   
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
            self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Splice Variants ###
            if self.opt['SpliceVar'] and splicevar:
                ## Setup search ##
                splicevar.sort()    #X#                splicevar = rje.sortUnique(splicevar)
                if self.list['Extract']:
                    self.deBug(self.list['Extract'])
                    self.list['Extract'] += splicevar
                    self.deBug(self.list['Extract'])
                self.log.printLog('#VAR','Looking for %s potential splice variants.' % (rje.integerString(len(splicevar))))
                start_pos = true_start
                for acc in splicevar[0:]:
                    ## Search index ##
                    ipos = rje.posFromIndex(acc,INDEX,start_pos,end_pos,re_index)   #X#,sortunique=True)
                    if ipos < 0:    # Could not find target!
                        missing.append(acc)
                        continue
                    ## Update ##
                    (line,start_pos) = rje.fileLineFromSeek(INDEX,ipos,reseek=False,next=False)
                    (matchacc,key,pos) = rje.matchExp(re_index,line)        #INDEX.readline())
                    if dat_dict[key].has_key(pos):  # Dictionary of (key:{dictionary of {positions:list of accs/ids}})
                        dat_dict[key][pos].append(acc)
                    else:
                        dat_dict[key][pos] = [acc]
                    if acc_dict.has_key(acc):   # Dictionary of (acc/id:{list of {key:pos}}) (single acc can have multiple entries)
                        acc_dict[acc].append({key:pos})  
                    else:                    
                        acc_dict[acc] = [{key:pos}]   
                    self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))),log=False,newline=False)
                self.log.printLog('\r#INDEX','Found index entries for %s of %s AccNum/ID. %s missing.' % (rje.integerString(len(acc_dict)),rje.integerString(len(acc_list)),rje.integerString(len(missing))))

            ### Missing acc have no entry in returned dictionary ###
            #X#for acc in missing:            
            #X#    self.log.printLog('#ACC','AccNum/ID "%s" missing from %s' % (acc,indexfile))
            
            ### Extract From UniProt using dictionaries ###
            _stage = 'Extract using Dictionaries'
            extract_dict = {}   # {{key:pos}:'ID (AccNum)'} = For matching input acclist to extracted entries
            for key in rje.sortKeys(dat_keys):
                unifile = dat_keys[key]
                UNIFILE = open(unifile,'r')
                ex = 0
                for pos in dat_dict[key]:
                    ipos = string.atol(pos)
                    # Get to correct position
                    rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[1]
                    # Read
                    if self._readSingleEntry(UNIFILE,logft=False,cleardata=cleardata):
                        # Check for replaced AccNum/ID
                        wanted_acc = dat_dict[key][pos]
                        for acc in wanted_acc:
                            accentry[acc] = self.list['Entry'][-1]
                        # Update
                        ex += 1
                        self.list['Entry'] = self.list['Entry'][:-1]     # Delete 
                    else:
                        bummer = rje.chomp(rje.fileLineFromSeek(UNIFILE,ipos,reseek=True,next=False)[0])
                        self.log.errorLog('%s rejected by _readSingleEntry() but explicitly selected for extraction!' % bummer,printerror=False)
                    self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile),log=False,newline=False)
                UNIFILE.close()
                self.log.printLog('\r#ACC','%s entries extracted from %s.' % (rje.integerString(ex),unifile))

            return (accentry)
            
        except:
            self.log.errorLog('Program error during rje_uniprot.accEntry()',quitchoice=True)
#########################################################################################################################
## End of SECTION III: UniProt Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: UniProtEntry Class                                                                                      # 
#########################################################################################################################
class UniProtEntry(rje.RJE_Object):     
    '''
    UniProt Entry Class. Author: Rich Edwards (2005).

    Info:str
    - Name = UniProt ID of Entry
    - Type = Preliminary or Standard
    - FullText = Full Text of Entry
    
    Opt:boolean
    - CC2FT = Extra whole-length features added for TISSUE and LOCATION [False]
    - ClearData = Whether to clear unprocessed Entry data (True) or (False) retain in Entry & Sequence objects [True]
    - InvMask = Whether to invert the masking and only retain maskft features [False]    
    - TMConvert = Whether to convert TOPO_DOM features, using first description word as Type [False]

    Stat:numeric
    - Length = Length of Sequence as annotated

    List:list
    - CaseFT = List of Features to make upper case with rest of sequence lower case []
    - Feature = List of feature dictionaries: [Type,Start,End,Desc]
    - MaskFT = List of Features to mask out []
    - PubMed = List of PubMed IDs (as strings)
    - Keywords = List of UniProt Keywords
    - Tissues = List of UniProt Tissues
    - Synonyms = List of Gene synonyms
    
    Dict:dictionary
    - Data = Dictionary of lists of UniProt data (Keys are line headers ID/AC/CC etc.)
    - DB = Specific extractions from DR lines for use in other programs. {DB:[AccNum/ID]}
    - Comments = Dictionary of comments: {Type:List of Comments}
    - DBLinks = List of Database Link dictionaries {Dbase,List of Details} for dblinks output

    Obj:RJE_Objects
    - Sequence = rje_sequence.Sequence object
    '''
    ### Attributes
    def shortName(self):    return '%s (%s)' % (self.obj['Sequence'].info['ID'],self.obj['Sequence'].info['AccNum'])
    def seqi(self,ikey): return self.obj['Sequence'].info[ikey]
#########################################################################################################################
    def isSpecies(self,spec=None,speclist=[]):  ### Returns True if entry corresponds to listed species
        '''
        Returns True if entry corresponds to listed species (or species code).
        >> spec:str = single species that MUST be be right
        >> speclist:list = can match any species in list
        '''
        if spec: return spec in [self.obj['Sequence'].info['Species'],self.obj['Sequence'].info['SpecCode']]
        for sp in speclist:
            if sp in [self.obj['Sequence'].info['Species'],self.obj['Sequence'].info['SpecCode']]: return True
        return False
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','FullText']
        self.statlist = ['Length']
        self.optlist = ['CC2FT','InvMask','TMConvert','ClearData']
        self.listlist = ['Feature','PubMed','Keywords','Tissues','Synonyms','MaskFT','CaseFT']
        self.dictlist = ['Data','Comments','DBLinks','DB']
        self.objlist = ['Sequence']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.obj['Sequence'] = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list)
        self.info['FullText'] = ''
        self.list['Feature'] = []   # List of features = {'Type':str,'Start':int,'End':int,'Desc':str}
        self.list['MaskFT'] = []
        self.setOpt({'ClearData':True})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'opt',['CC2FT','InvMask','TMConvert','ClearData'])
                self._cmdReadList(cmd,'list',['CaseFT','MaskFT'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Attribute Processing
#########################################################################################################################
    def _uniParse(self,key):    ### Parses list of elements from self.dict['Data'] (and pops) 
        '''
        Parses list of elements from self.dict['Data'] (and pops).
        >> key:str = Key of UniProt entry type
        << List of matched elements or False if failure.
        '''
        try:### ~ [1] Check key and parse data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if key not in self.dict['Data'].keys(): return False
            return rje.matchExp(uniparse[key],self.dict['Data'][key][0])
            #!# This no longer pops elements from dictionary: check backwards compatibility #!#
        except:
            self.log.errorLog('UniProtEntry: Cataclysmic error during _uniParse(%s)!' % key)
            return False
#########################################################################################################################
    def process(self,logft=True,cleardata=None):  ### Extract Details from self.dict['Data'] to Sequence object
        '''
        Extract Details from self.dict['Data'] to Sequence object.
        >> logft:boolean = whether to write number of features to log file
        >> cleardata=Whether to clear self.dict['Data'] after processing (to save memory) [True]
        << True if OK, False if not.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            if cleardata == None: cleardata = self.opt['ClearData']
            seq_info = ['Name','Type','Description','Sequence','ID','AccNum','DBase','Gene','Species','SpecCode','Format','TaxaID']
            seqi = self.obj['Sequence'].info
            for i in seq_info:
                if i not in seqi: seqi[i] = ''
            seqi['Type'] = 'Protein'
            #X#print '\nAcc:', self.dict['Data']['AC']
            _stage = 'Compress Certain uniprot lists'
            for key in ['DE','GN','OC','AC']:
                if self.dict['Data'].has_key(key): self.dict['Data'][key] = [string.join(self.dict['Data'][key])]
                                    
            ### Basic Sequence Details ###
            _stage = 'Basic Details (ID)'
            parse = self._uniParse('ID')
            if parse:
                (self.info['Name'],self.info['Type'],self.stat['Length']) = parse[:3]
                self.stat['Length'] = string.atoi(self.stat['Length'])
            else: self.stat['Length']
            self.info['ID'] = seqi['ID'] = self.info['Name']

            _stage = 'AccNum (AC)'
            full_acc = string.split(string.join(self.dict['Data']['AC']))
            self.obj['Sequence'].list['Secondary ID'] = string.split(string.join(full_acc,''),';')[1:-1]
            parse = self._uniParse('AC')
            if parse: seqi['AccNum'] = parse[0]
            if string.split(self.info['Name'])[0] == seqi['AccNum']:
                seqi['ID'] = seqi['AccNum']
                seqi['DBase'] = 'custom'
            elif self.info['Name'].find(seqi['AccNum']) == 0: seqi['DBase'] = 'trembl'
            else: seqi['DBase'] = 'sprot'
            if seqi['ID'] != seqi['AccNum']: self.info['Name'] = '%s__%s' % (seqi['ID'],seqi['AccNum'])

            _stage = 'Description (DE)'
            parse = self._uniParse('DE')
            if parse: seqi['Description'] = parse[0]
                
            _stage = 'Gene (GN)'
            parse = self._uniParse('GN')
            if parse: seqi['Gene'] = parse[0]
            if self.dict['Data'].has_key('GN'):
                if rje.matchExp(uniparse['SY'],self.dict['Data']['GN'][0]):
                    syn = string.split(rje.matchExp(uniparse['SY'],self.dict['Data']['GN'][0])[0],';')[0]
                    self.list['Synonyms'] = string.split(syn,', ')
                    
            _stage = 'Species (OS)'
            parse = self._uniParse('OS')
            if parse: seqi['Species'] = parse[0]
            parse = self._uniParse('OX')
            if parse: seqi['TaxaID'] = parse[0]
            seqi['SpecCode'] = seqi['ID'][seqi['ID'].find('_')+1:]
            
            _stage = 'Name'
            if seqi['DBase'] == 'trembl' and seqi['Gene'] != 'None':
                for g in seqi['Gene'][0:]:
                    if not rje.matchExp('([A-Za-z0-9_-])',g) and g not in ['.','#']: seqi['Gene'] = string.replace(seqi['Gene'],g,'')
                seqi['ID'] = '%s_%s' % (seqi['Gene'].lower(),seqi['SpecCode'])
            if seqi['ID'] == seqi['AccNum']:
                seqi['Name'] = '%s %s' % (seqi['ID'],seqi['Description'])
                seqi['Format'] = 'gnspec'
            else:
                seqi['Name'] = '%s__%s %s' % (seqi['ID'],seqi['AccNum'],seqi['Description'])
                seqi['Format'] = 'gn_sp__acc'
            
            _stage = 'Sequence'
            if self.dict['Data'].has_key('SEQ'): seqi['Sequence'] = self.dict['Data'].pop('SEQ')[0]
            seqi['Sequence'] = re.sub('\s+','',seqi['Sequence']).upper()
            if not self.stat['Length']: self.stat['Length'] = len(seqi['Sequence'])

            ### Features ###
            _stage = 'Features'
            if self.dict['Data'].has_key('FT'):
                for ft in self.dict['Data']['FT']:
                    # Remove '?'
                    while rje.matchExp('(\?\d)',ft) or rje.matchExp('(\d\?)',ft): ft = string.replace(ft,'?','')
                    ft = string.replace(ft,'?','0')
                    parse = rje.matchExp(uniparse['FT'],ft)
                    parse_nodesc = rje.matchExp('(\S+)\s+<*(\d+)\s+>*(\d+)\.*',ft)
                    parse_onepos = rje.matchExp('(\S+)\s+(\d+)\.*\s+(\S+)',ft)
                    parse_cntd = rje.matchExp('\s+(\S.*)$',ft)
                    if parse:
                        ftdic = {
                            'Type' : parse[0],
                            'Start' : string.atoi(parse[1]),
                            'End' : string.atoi(parse[2]),
                            'Desc' : parse[3]
                            }
                        if rje.matchExp(string.join(['(\S+)','<(\d+)','>*(\d+)\.*','(\S.+)\s*$'], '\s+'),ft) or rje.matchExp(string.join(['(\S+)','<*(\d+)','>(\d+)','(\S.+)\s*$'], '\s+'),ft):
                            ftdic['Desc'] = ftdic['Desc'] + ' (Truncated?)'
                        ftdic['Desc'] = re.sub('\s+',' ',ftdic['Desc'])
                        if self.opt['TMConvert'] and ftdic['Type'] == 'TOPO_DOM':
                            ftdic['Type'] = string.strip(string.split(ftdic['Desc'])[0].upper(),'.')
                            ftdic['Desc'] = 'TOPO_DOM %s' % ftdic['Desc']
                        self.list['Feature'].append(ftdic)
                    elif parse_nodesc:
                        parse = parse_nodesc
                        if parse:
                            ftdic = {
                                'Type' : parse[0],
                                'Start' : string.atoi(parse[1]),
                                'End' : string.atoi(parse[2]),
                                'Desc' : parse[0]
                                }
                            if rje.matchExp('(\S+)\s+<(\d+)\s+>*(\d+)\.*',ft) or rje.matchExp('(\S+)\s+<*(\d+)\s+>(\d+)\.*',ft):
                                ftdic['Desc'] = ftdic['Desc'] + ' (Truncated?)'
                            self.list['Feature'].append(ftdic)
                    elif parse_onepos:
                        parse = parse_onepos
                        if parse:
                            ftdic = {
                                'Type' : parse[0],
                                'Start' : string.atoi(parse[1]),
                                'End' : string.atoi(parse[1]),
                                'Desc' : parse[2]
                                }
                            self.list['Feature'].append(ftdic)
                    elif parse_cntd:
                        try:
                            ftdic = self.list['Feature'][-1]
                            ftdic['Desc'] = re.sub('\s+',' ','%s %s' % (ftdic['Desc'],parse_cntd[0]))
                        except: self.printLog('#ERR','No feature parsed for addition of "%s".' % ft)
                    else: self.log.printLog('#ERR','Cannot parse feature details from %s.' % ft)

            ### Tissues (RC) ###
            _stage = 'Tissues'
            if self.dict['Data'].has_key('RC'):
                tissues = []
                self.list['Tissues'] = []
                for rc in self.dict['Data']['RC']:
                    parse = rje.matchExp(uniparse['RC'],rc)
                    if parse: tissues += string.split(parse[0],', ')
                for tissue in tissues:
                    if tissue[:4] == 'and ': self.list['Tissues'].append(tissue[4:])
                    else: self.list['Tissues'].append(tissue)

            ### Keywords (KW) ###
            _stage = 'Keywords'
            if self.dict['Data'].has_key('KW'):
                keywords = string.join(self.dict['Data']['KW'])
                if keywords[-1:] == '.': keywords = keywords[:-1]    # Remove full stop
                self.list['Keywords'] = string.split(keywords,'; ')

            ### References (RX) ###
            _stage = 'PubMed IDs'
            if self.dict['Data'].has_key('RX'):
                self.list['PubMed'] = []
                for rx in self.dict['Data']['RX']:
                    parse = rje.matchExp(uniparse['RX'],rx)
                    if parse: self.list['PubMed'].append(parse[0])

            ### Comments (CC) ###
            _stage = 'Comments'
            if self.dict['Data'].has_key('CC'):
                self.dict['Comments'] = {}  # Dictionary of comments: {Type:List of Comments}
                for cc in self.dict['Data']['CC']:
                    if cc.find('-----') == 0: break
                    csplit = string.split(cc[4:],': ')
                    ctype = csplit[0]
                    cdetail = string.join(csplit[1:],': ')
                    if self.dict['Comments'].has_key(ctype): self.dict['Comments'][ctype].append(cdetail)
                    else: self.dict['Comments'][ctype] = [cdetail]

            ### Database Links ### (DR)
            _stage = 'Database Links'
            if self.dict['Data'].has_key('DR'):
                self.dict['DBLinks'] = {}   # Database Link dictionary {Dbase,List of Details}
                for dr in self.dict['Data']['DR']:
                    if rje.matchExp(uniparse['DR'],dr):
                        (ctype,cdetail) = rje.matchExp(uniparse['DR'],dr)
                        if self.dict['DBLinks'].has_key(ctype): self.dict['DBLinks'][ctype].append(cdetail)
                        else: self.dict['DBLinks'][ctype] = [cdetail]
                        self.specialDB(ctype,cdetail)   # Extracts specific information to self.dict['DB']
                    elif rje.matchExp('^(Entrez Gene); (\S+)$',dr):
                        (ctype,cdetail) = rje.matchExp('^(Entrez Gene); (\S+)$',dr)
                        self.specialDB(ctype,cdetail)   # Extracts specific information to self.dict['DB']

            ### CC to FT ###
            _stage = 'End'
            if self.opt['CC2FT']: self.cc2ft()

            ### FT Masking/Case Change ###
            if self.list['MaskFT']: self.maskFT(self.list['MaskFT'],inverse=self.opt['InvMask'])
            if self.list['CaseFT']: self.caseFT(self.list['CaseFT'])
            #X#for k in ['Synonyms','PubMed','Keywords','Tissues']:
            #X#    print k, self.list[k]
            #X#for k in ['Comments','DBLinks']:
            #X#    print k, self.dict[k]         

            ### Cleanup ###
            #X#for key in self.dict['Data'].keys():
            #X#    if key not in useful_data:
            #X#        self.dict['Data'].pop(key)      # Save memory!
            #X#print self.dict['Data']
            #X#print self.info, self.obj['Sequence'].info
            if cleardata: self.dict['Data'] = {}    # Save memory!
            else: self.obj['Sequence'].dict['UniDAT'] = self.dict['Data']
            if logft:
                self.printLog('#FT','%d features for %s.' % (len(self.list['Feature']),self.obj['Sequence'].info['AccNum']))
            return True
        except:
            self.log.errorLog('Cataclysmic error during UniProtEntry: process() %s!' % _stage)
            return False
#########################################################################################################################
    def specialDB(self,dbase,details):  ### Extracts specific information to self.dict['DB']
        '''
        Extracts specific information to self.dict['DB'].
        >> dbase:str = Database identifier extracted from DR line of DAT file - '^(\S+);\s+(\S.+)$'[0]
        >> details:str = Database links extracted from DR line of DAT file - '^(\S+);\s+(\S.+)$'[1]
        '''
        try:    ### Each Database is dealt with individually. If the database is not here, nothing will happen ###
            ###~UniProt (From IPI DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   UniProtKB/Swiss-Prot; O95793-1; STAU1_HUMAN; M.
            if dbase == "UniProtKB/Swiss-Prot":
                try: (uni_sv,uni_id) = rje.matchExp('^(\S+); (\S+);',details)
                except: return
                uni_acc = string.split(uni_sv,'-')[0]
                for db in ['UniAccNum','UniID','UniSV']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                self.dict['DB']['UniAccNum'].append(uni_acc)
                self.dict['DB']['UniID'].append(uni_id)
                if uni_acc.find('-') > 0: self.dict['DB']['UniSV'].append(uni_sv)

            ###~EnsEMBL (From IPI/UniProt DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   ENSEMBL; ENSP00000346163; ENSG00000124214; -.       = IPI
            ## DR   Ensembl; ENSG00000112530; Homo sapiens.             = UniProt
            if dbase.lower() == 'ensembl':
                for db in ['EnsG','EnsP']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                if rje.matchExp('(ENSP\S+);',details): self.dict['DB']['EnsP'].append(rje.matchExp('(ENSP\S+);',details)[0])
                if rje.matchExp('(ENSG\S+);',details): self.dict['DB']['EnsG'].append(rje.matchExp('(ENSG\S+);',details)[0])

            ###~RefSeq/NCBI (From IPI DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   REFSEQ_REVIEWED; NP_059347; GI:82659087; -.
            ## DR   REFSEQ_PROVISIONAL; NP_004520; GI:4758720; -.
            ## DR   REFSEQ_VALIDATED; NP_001014313; GI:62122851; -.
            if dbase.find('REFSEQ') == 0:
                if not self.dict['DB'].has_key('RefSeq'): self.dict['DB']['RefSeq'] = []
                self.dict['DB']['RefSeq'].append(rje.matchExp('(\S+_\S+);',details)[0])

            ###~HGNC (From IPI/UniProt DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   HGNC; 11370; STAU1; -.      = IPI
            ## DR   HGNC; HGNC:19152; PACRG.    = UniProt
            if dbase == 'HGNC':
                for db in ['Symbol','HGNC']:
                    if not self.dict['DB'].has_key(db): self.dict['DB'][db] = []
                if rje.matchExp('^HGNC:(\d+); (\S+).',details): (hgnc,symbol) = rje.matchExp('^HGNC:(\d+); (\S+).',details)
                else: (hgnc,symbol) = rje.matchExp('^(\d+); (\S+);',details)
                self.dict['DB']['HGNC'].append(hgnc)
                self.dict['DB']['Symbol'].append(symbol)

            ###~Entrez Gene (From IPI DAT)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## DR   Entrez Gene; 6780; STAU1; -.
            if dbase == 'Entrez Gene':
                if not self.dict['DB'].has_key('Entrez'): self.dict['DB']['Entrez'] = []
                self.dict['DB']['Entrez'].append(rje.matchExp('^(\d+);',details)[0])

        except: self.errorLog('Problem with UniProtEntry.specialDB(%s)' % dbase)
#########################################################################################################################
    def cc2ft(self):    ### Adds extra full-length features based on LOCATION and TISSUE
        '''Adds extra full-length features based on LOCATION and TISSUE.'''
        try:
            ### Setup ###
            pos = (1,self.obj['Sequence'].aaLen())
            go = []
            if self.dict['DBLinks'].has_key('GO'):
                go = self.dict['DBLinks']['GO'][0:]
                for g in go[0:]:
                    if g.find('; C:') < 0:
                        go.remove(g)
            cc = 'SUBCELLULAR LOCATION'
            if self.dict['Comments'].has_key(cc):
                go = string.split(string.join(self.dict['Comments'][cc],'; '),'; ') + go
            ftlist = {'TISSUE':self.list['Tissues'],'LOCATION':go}
            ### Add FT ##
            for type in ftlist.keys():
                for desc in ftlist[type]:
                    if {'Type':type,'Desc':desc,'Start':pos[0],'End':pos[1]} not in self.list['Feature']:
                        self.list['Feature'].append({'Type':type,'Desc':desc,'Start':pos[0],'End':pos[1]})
            return True                
        except:
            self.log.errorLog('UniProtEntry.cc2ft() has gone wrong!')
            return False
#########################################################################################################################
    def maskFT(self,types=['EM'],inverse=False,mask='X',log=True):   ### Masks given feature types
        '''
        Masks given feature types.
        >> types:list of str [['EM']] = types of feature to mask
        >> inverse:bool [False] = whether to mask all sequence *except* listed types
        >> mask:str ['X'] = character to use for masking
        >> log:bool [True] = whether to log the affects of masking
        '''
        try:
            ### Setup ###
            oldseq = self.obj['Sequence'].info['Sequence'][0:]
            newseq = mask * len(oldseq)
            mx = 0
            prex = oldseq.count(mask)
            if inverse:
                (newseq,oldseq) = (oldseq,newseq)
            ### Mask ###
            for ft in self.list['Feature']:
                if ft['Type'] in types:
                    oldseq = oldseq[:ft['Start']-1] + newseq[ft['Start']-1:ft['End']] + oldseq[ft['End']:]
                    mx += 1
            ### Update ###
            self.obj['Sequence'].info['Sequence'] = oldseq
            maskx = oldseq.count(mask) - prex
            if log:
                mtxt = self.info['Name']
                if inverse:
                    mtxt += ' inverse'
                self.log.printLog('#MASK','%s masked %d features. (%d %s added.)' % (mtxt,mx,maskx,mask))
        except:
            self.log.errorLog('Problem masking %s features from %s' % (types,self.shortName()))
#########################################################################################################################
    def caseFT(self,types=[]):   ### Masks given feature types
        '''
        Masks given feature types.
        >> types:list of str [] = types of feature to be upper case
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sequence = self.obj['Sequence'].info['Sequence'][0:].lower()
            ### ~ [2] ~ Change Case ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for ft in self.list['Feature']:
                if ft['Type'] in types:
                    sequence = sequence[:ft['Start']-1] + sequence[ft['Start']-1:ft['End']].upper() + sequence[ft['End']:]
            self.obj['Sequence'].dict['Case'] = rje_sequence.caseDict(sequence)
        except:
            self.log.errorLog('Problem masking %s features from %s' % (types,self.shortName()))
#########################################################################################################################
    def orderFT(self):  ### Orders features by start, end, type
        '''Orders features by start, end, type.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newft = []
            ### ~ [2] ~ Order ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for ft in self.list['Feature']:
                i = 0
                while i < len(newft):
                    if ft['Start'] < newft[i]['Start']: break
                    elif ft['Start'] > newft[i]['Start']: i += 1; continue
                    if ft['End'] < newft[i]['End']: break
                    elif ft['End'] > newft[i]['End']: i += 1; continue
                    if ft['Type'] < newft[i]['Type']: break
                    elif ft['Type'] > newft[i]['Type']: i += 1; continue
                    i += 1
                newft.insert(i,ft)
            self.list['Feature'] = newft
        except: self.errorLog('Problem ordering %s features' % self.info['Name'])            
#########################################################################################################################
    ### <3> ### UniProt Conversion and Saving                                                                           #
#########################################################################################################################
    def uniProtFromSeq(self,seq=None,sequence='',name='',data={},ft=[]): ### Converts into UniProtEntry object (self!)
        '''
        Converts into UniProtEntry object (self!).
        >> seq:rje_sequence.Sequence object [None]
        >> sequence:str = alternative sequence data (will be converted to Sequence object!) ['']
        >> name:str = alternative sequence name (will be converted to Sequence object!) ['']
        >> data:dict = dictionary of UniProt data with {keys ID/AC/OS etc: [list of lines]} [{}]
        >> ft:list = list of ftdic dictionaries of features {'Type/Desc':str,'Start/End':int} [[]]
        << returns self if successful or None if fails
        '''
        try:
            ### Setup ###
            if not seq and (sequence and name):
                seq = rje_sequence.Sequence(log=self.log)
                seq.info['Name'] = name
                seq.info['Sequence'] = sequence.upper()  #!# Change at some point to allow mixed case!
                seq.info['Type'] = 'Protein'
                seq.extractDetails()    #!#gnspacc=self.opt['GeneSpAcc'])
            if not seq:
                seq = self.obj['Sequence']
            if not seq:
                raise ValueError, 'No sequence information given'
            self.obj['Sequence'] = seq

            ### Update self ###
            for key in data:
                self.dict['Data'][key] = data[key]
            self.list['Feature'] += ft
            if seq.info['DBase'] != 'trembl':
                self.dict['Data']['ID'] = ['%s     %s;   %d AA.\n' % (seq.info['ID'],seq.info['Type'],seq.aaLen())]
            elif seq.info['SpecCode'] not in ['None','UNK']:
                self.dict['Data']['ID'] = ['%s_%s     %s;   %d AA.\n' % (seq.info['AccNum'],seq.info['SpecCode'],seq.info['Type'],seq.aaLen())]
            else:
                self.dict['Data']['ID'] = ['%s     %s;   %d AA.\n' % (seq.info['AccNum'],seq.info['Type'],seq.aaLen())]
            if self.dict['Data'].has_key('AC'):
                self.dict['Data']['AC'] = ['%s;' % seq.info['AccNum']] + self.dict['Data']['AC']
            else:
                self.dict['Data']['AC'] = ['%s;' % seq.info['AccNum']]
            if seq.info['Description'].lower() != 'none':
                self.dict['Data']['DE'] = [seq.info['Description']]
            else:
                self.dict['Data']['DE'] = ['']
            if seq.info['Species'] not in ['None','Unknown',seq.info['SpecCode'],'']:
                if self.dict['Data'].has_key('OS'):
                    self.dict['Data']['OS'] = ['%s.' % seq.info['Species']] + self.dict['Data']['OS']
                else:
                    self.dict['Data']['OS'] = ['%s.' % seq.info['Species']]
            dt = string.split(time.ctime())
            if self.dict['Data'].has_key('DT'):
                self.dict['Data']['DT'] = ['%s-%s-%s, generated by rje_uniprot' % (dt[2],dt[1].upper(),dt[-1])] + self.dict['Data']['DT']
            else:
                self.dict['Data']['DT'] = ['%s-%s-%s, generated by rje_uniprot' % (dt[2],dt[1].upper(),dt[-1])]
            if self.dict['Data'].has_key('CC'):
                self.dict['Data']['CC'] = ['-!- Entry generated by rje_uniprot %s' % time.ctime()] + self.dict['Data']['CC']
            else:
                self.dict['Data']['CC'] = ['-!- Entry generated by rje_uniprot %s' % time.ctime()]

            #X#self.deBug(self.dict['Data'])
            if self.process(logft=False,cleardata=False):
                return self
            return None
        except:
            self.log.errorLog('UniProtEntry.uniProtFromSeq() has gone wrong.',quitchoice=True)
            return None
#########################################################################################################################
## End of SECTION III : UniProtEntry Class                                                                              #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: General UniProt Methods                                                                                 #
#########################################################################################################################
def downloadUniProt(callobj):   ### Downloads the UniProt database using the attributes of callobj
    '''Downloads the UniProt database using the attributes of callobj.'''
    try: ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        outdir = callobj.info['UniPath']
        mydir = os.path.abspath(os.curdir)
        if not os.path.exists(outdir): rje.mkDir(callobj,outdir)
        ftproot = 'ftp://ftp.ebi.ac.uk/pub/databases/uniprot/current_release/'
        dbtxt = 'UniProt from %s -> %s' % (ftproot,outdir)
        if callobj.stat['Interactive'] >= 0 and not rje.yesNo('Download %s?' % dbtxt): return False
        elements = ['knowledgebase/complete/*dat.gz','README','relnotes.txt']

        ### ~ [2] Download files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        callobj.log.printLog('#DBASE','Downloading %s: %d elements.' % (dbtxt,len(elements)))
        for element in elements:
            callobj.log.printLog('#DBASE','Downloading %s' % (element))
            os.chdir(outdir)
            os.system('wget %s' % rje.makePath(ftproot+element,wholepath=True))
            os.chdir(mydir)
        return True
    except: callobj.log.errorLog('Major error during UniProt download')
#########################################################################################################################
def processUniProt(callobj,makeindex=True,makespec=True,makefas=True):     ### Processes UniProt making index file and spectable as appropriate
    '''Processes UniProt making index file and spectable as appropriate.'''
    try:### ~ [1] Setup Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not makeindex and not makespec and not makefas:
            return callobj.log.printLog('#DB','No call to make UniProt index/spectable/fasta files.')
        unipath = callobj.info['UniPath']
        try: indexfile = unipath + callobj.info['DBIndex']
        except: indexfile = unipath + 'uniprot.index'
        datfiles = glob.glob('%s*.dat' % unipath)
        datfiles.sort()
        if not datfiles: return callobj.log.errorLog('No *.dat files found in %s!' % unipath,printerror=False)
        ## ~ [1a] Check index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makeindex or callobj.stat['Interactive'] < 1 or rje.yesNo('Check/Make UniProt index file?'):
            if not os.path.exists(indexfile) or callobj.opt['Force']: makeindex = True
            else:
                makeindex = False
                for dat in datfiles:
                    if rje.isYounger(indexfile,dat) != indexfile:
                        makeindex = True
                        break
        ## ~ [1b] SpecTable Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        # - can then use "grep 'Metazoa' uniprot.spec.tdt | gawk '{print $1 }' | sort -u" to extract, e.g. Metazoan species codes
        spectable = unipath + 'uniprot.spec.tdt'
        spectemp = unipath + 'uniprot.spec.txt'
        if makespec or (callobj.stat['Interactive'] > 0 and rje.yesNo('Check/Make UniProt species file?')):
            if callobj.opt['Force'] or not os.path.exists(spectable): makespec = True
            else:
                makespec = False
                for dat in datfiles:
                    if rje.isYounger(spectable,dat) != spectable:
                        makespec = True
                        break
        ## ~ [1c] Fasta Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makefas or (callobj.stat['Interactive'] > 0 and rje.yesNo('Reformat UniProt to fasta files?')):
            if callobj.opt['Force']: makefas = True
            else:
                makefas = False
                for dat in datfiles:
                    fas = rje.baseFile(dat) + '_all.fas'
                    if not os.path.exists(fas) or rje.isYounger(fas,dat) != fas:
                        makefas = True
                        break
        ## ~ [1d] Abort if no need to process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not makeindex and not makespec and not makefas:
            return callobj.log.printLog('#DB','No need to make UniProt index/spectable/fasta files. (Force=F)')
        callobj.deBug('Index: %s; Spec: %s; Fas: %s' % (makeindex,makespec,makefas))

        ### ~ [2] Setup Output and Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ix = -1 # No. index lines
        tx = 0  # No. taxa lines
        sx = 0  # No. seqs
        dx = 0  # No. dat files
        if makeindex:
            INDEX = open('%s.head' % indexfile,'w')
            for d in range(len(datfiles)):
                dat = datfiles[d]
                INDEX.write('#%d=%s\n' % (d+1,os.path.basename(dat)))
            INDEX.close()
            INDEX = open('%s.temp' % indexfile,'w')     # Clear it
        if makespec: SPEC = open(spectemp,'w')          # Clear it
        if makefas:
            (append,callobj.opt['Append']) = (callobj.opt['Append'],False)
            for dat in datfiles: rje.backup(callobj,rje.baseFile(dat) + '_all.fas')
            callobj.opt['Append'] = append

        ### ~ [3] Process Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        wanted = ['ID','AC','//']
        if makespec: wanted += ['OS','OC','OX']
        if makefas: wanted += ['DE','  ']
        for dat in datfiles:
            ## ~ [3a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fasfile = rje.baseFile(dat) + '_all.fas'
            DAT = open(dat, 'r')
            ix = seq_index = -1
            dx += 1
            dbtext = 'Processing %d of %d files (%s):' % (dx,len(datfiles),dat)
            id = acc = code = species = seq = taxid = desc = ''
            taxonomy = []
            file_pos = 0
            ## ~ [3b] Read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while 'Reading File':
                line = DAT.readline()
                if not line:
                    DAT.close()
                    break
                ix += 1
                if line[:2] not in wanted: continue
                ## ~ [3c] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line[0:2] == 'ID':
                    id = string.split(line)[1]
                    try: code = string.split(id,'_')[1]
                    except: pass
                    seq_index = file_pos
                elif line[0:2] == 'AC':
                    extras = rje.matchExp('^AC\s+(\S.+)\s*$',line)[0]
                    if acc: acc = '%s %s' % (acc,extras)
                    else: acc = extras
                elif line[:2] == 'DE':
                    if desc: desc += rje.chomp(line[4:])
                    else: desc = rje.chomp(line[5:])
                elif line[0:2] == 'OS' and not species:
                    try: species = rje.matchExp('^OS\s+(\S[A-Za-z0-9\-\s:\'\/#]+\S)\s+\(+',line)[0]
                    except:
                        try: species = rje.matchExp('^OS\s+(\S[A-Za-z0-9\-\s:\'\/#]+\S)[\.,]',line)[0]
                        except:
                            try: species = rje.matchExp('^OS\s+(\S.+irus)',line)[0]
                            except: species = rje.matchExp('^OS\s+(\S.+\S)',line)[0]
                    if not code: code = rje_sequence.getSpecCode(species)
                elif line[0:2] == 'OC':
                    taxonomy = taxonomy + string.split(re.sub('\s+','',rje.matchExp('^OC\s+(\S.+)$',line)[0]),';')
                    if '' in taxonomy: taxonomy.remove('')
                elif line[0:2] == 'OX':
                    taxid = rje.matchExp('^OX\s+NCBI_TaxID=(\d+)',line)[0]
                    taxonomy.append(species)
                elif line[:2] == '  ': seq += string.replace(rje.chomp(line[5:]),' ','')
                elif line[0:2] == '//':
                    sx += 1
                    file_pos = DAT.tell()
                    if makefas and seq: open(fasfile,'a').write('>%s__%s %s\n%s\n' % (id,string.split(acc,';')[0],desc,seq))
                    if code and makespec:
                        tx += 1
                        SPEC.write('%s\t:%s:\t:%s:\n' % (code, string.join(taxonomy,':'), taxid))
                    if seq_index >= 0 and makeindex:
                        ilist = string.split(string.join(string.split(acc) + [id],''),';')
                        while '' in ilist: ilist.remove('')
                        for acc in ilist: INDEX.write('%s;%d:%d\n' % (acc,dx,seq_index))
                    id = acc = code = species = seq = taxid = desc = ''
                    taxonomy = []
                    seq_index = -1
                    callobj.log.printLog('\r#DB','%s %s lines; %s entries' % (dbtext,rje.integerString(ix),rje.integerString(sx)),log=False,newline=False)
            callobj.log.printLog('\r#DB','%s %s lines; %s entries' % (dbtext,rje.integerString(ix),rje.integerString(sx)))                    
    
        ### ~ [4] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        dbtext = 'Processed %d *.dat files:' % len(datfiles)
        callobj.log.printLog('#DB','%s %s lines; %s entries.' % (dbtext,rje.integerString(ix),rje.integerString(sx)))
        ### ~ [4a] Sort out index files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makeindex:
            INDEX.close()

            ### Subsort within subsets ##
            subset = 3  # Cannot use 2 even though I want to because of fucking dumb UNIX sort P_ problem. Grrrr!
            TMP = open('%s.temp' % indexfile,'r')
            subindex = []   # List of temp subindex files
            ## Read first sort and split into subsorts ##
            tmplines = []
            reading = True
            px = 0
            while reading:
                line = TMP.readline()
                if line in ['',None]: break
                px += 1
                start = 'index_%s.tmp' % line[:subset]
                if start not in subindex: subindex.append(start)
                open(start,'a').write(line)
                callobj.log.printLog('\r#INDEX','Processing index lines: %s' % rje.integerString(px),log=False,newline=False)
            callobj.log.printLog('\r#INDEX','Processed %s index lines into sub-files' % rje.integerString(px),log=False,newline=False)
            TMP.close()

            ### Process SubFiles ###
            subindex.sort()
            tmplines = callobj.loadFromFile('%s.head' % indexfile)
            INDEX = open(indexfile,'w')
            INDEX.write(string.join(tmplines,''))
            for file in subindex:
                TMP = open(file,'r')
                tmplines = []
                reading = True
                while reading:
                    line = TMP.readline()
                    if line in ['',None]: break
                    tmplines.append(line)
                indexdict = {}
                pxx = len(tmplines)
                px = 0
                while tmplines:
                    px += 50.0
                    indata = string.split(rje.chomp(tmplines.pop(0)),';')
                    if indexdict.has_key(indata[0]): indexdict[indata[0]] 
                    else: indexdict[indata[0]] = []
                    for idat in indata[1:]:
                        if idat not in indexdict[indata[0]]: indexdict[indata[0]].append(idat)
                    callobj.log.printLog('\r#INDEX','Processing %s index lines: %.2f%%' % (file,(px/pxx)),log=False,newline=False)
                for key in rje.sortKeys(indexdict):
                    px += 50.0
                    outdata = [key,string.join(indexdict.pop(key),';')]
                    INDEX.write('%s\n' % string.join(outdata,';'))
                    callobj.log.printLog('\r#INDEX','Processing %s index lines: %.2f%%' % (file,(px/pxx)),log=False,newline=False)
                callobj.log.printLog('\r#INDEX','Processing %s index lines: Complete.' % file,log=False)
                TMP.close()
                os.unlink(file)
            INDEX.close()
            os.unlink('%s.temp' % indexfile)
            os.unlink('%s.head' % indexfile)
            callobj.log.printLog('\r#INDEX','Generation of UniProt index complete.')
                        
            ### Old way using UNIX sort, which is quicker but doesn't work!! ###
            if 'you_want_to_do_it_the_quick_but_dangerous_way_with_UNIX_sort' in callobj.cmd_list:               
                if callobj.opt['Win32']: callobj.log.printLog('#DB','Cannot cleanup %s.temp with UNIX sort! Talk to Rich!' % indexfile)
                else:
                    callobj.log.printLog('#DB','Cleanup of %s.temp with UNIX sort ...' % indexfile,log=False)
                    if not os.system('sort %s.temp > %s.sort' % (indexfile,indexfile)):
                        os.unlink('%s.temp' % indexfile)
                        callobj.log.printLog('#DB','%s.temp sorted with UNIX sort. Concatenating...' % indexfile,log=False)
                        if not os.system('cat %s.head %s.sort > %s' % (indexfile,indexfile,indexfile)):
                            os.unlink('%s.sort' % indexfile)
                            os.unlink('%s.head' % indexfile)
                            callobj.log.printLog('#DB','Sort and Concatenation complete. Temporary files deleted. %s created.' % indexfile)
                        else: callobj.log.printLog('#ERR','Concatenation failed! Temporary files not deleted. %s not created.' % indexfile)                           
                    else: callobj.log.printLog('#ERR','Sorting of %s.temp with UNIX sort Failed!' % indexfile)

        ### ~ [4a] Sort out spectable files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if makespec:
            SPEC.close()
            if callobj.opt['Win32']:
                callobj.log.printLog('#DB','Cannot cleanup %s with UNIX sort. Very large file!' % spectemp)
                callobj.log.printLog('#DB','Use: "sort -u %s | grep -v \'^9\' > %s"' % (spectemp,spectable))
            else:
                callobj.log.printLog('#DB','Cleanup of %s with UNIX sort -u > %s ...' % (spectemp,spectable),log=False,newline=False)
                if not os.system('sort -u %s | grep -v \'^9\' > %s' % (spectemp,spectable)):
                    callobj.log.printLog('\r#DB','Cleanup of %s with UNIX sort -u > %s complete.' % (spectemp,spectable))
                    if callobj.stat['Interactive'] < 0 or rje.yesNo('Delete %s?' % spectemp): os.unlink(spectemp)
                else:
                    callobj.log.printLog('#DB','Cleanup of %s with UNIX sort failed. Very large file!' % spectemp)
                    callobj.log.printLog('#DB','Tried: "sort -u %s | grep -v \'^9\' > %s"' % (spectemp,spectable))
        
    except: callobj.log.errorLog('Error in processUniProt()',printerror=True)
#########################################################################################################################
## End of SECTION IV : Generic Module Methods                                                                           #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION V: MAIN PROGRAM                                                                                             #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try:        
        uniprot = UniProt(mainlog,cmd_list)
        uniprot.run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION V                                                                                                    #
#########################################################################################################################

