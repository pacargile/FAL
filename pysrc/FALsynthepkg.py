
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os,sys,glob,time,filecmp
import subprocess
import cStringIO
import shutil
from scipy import io as spIO
import numpy as np

from FALselMOD import selmod
import FALGlue

import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)

__all__ = ["synthepkg"]


class synthe(object):
    def __init__(self,ID=None,verbose=False,clobber=False,starpars=None):
        # define a unique job ID if not already defined
        # must be equal to or less than 8-digit int/string.
        if ID != None:
            self.ID = str(ID).rjust(8,str(0))
        else:
            import uuid
            self.ID = uuid.uuid4().hex[:8]
        with open('ID.asc','w') as f:
            f.write(self.ID)

        # if clobber=True, then remove all old local and
        # memory files
        if clobber == True:
            # first delete any *_* files in memory
            memfiles = glob.glob('/dev/shm/FAL/{0}/fort.*'.format(self.ID))
            for mf in memfiles:
                os.unlink(mf)
            # clear local directory of broken links
            locfiles = glob.glob('./*')
            for lf in locfiles:
                try:
                    os.stat(lf)
                except OSError:
                    os.unlink(lf)

        # define if writes to /dev/null or StdOUT
        if verbose == True:
            self.FNULL = sys.stdout
        elif verbose == False:
            self.FNULL = open(os.devnull, 'w')
        else:
            raise ValueError('Do not understand verbose flag')

        if starpars == None:
            raise ValueError('Must provide a starpars dictionary')
        else:
            self.starpars = starpars

        # Set up exec directory and data directory
        # self.HOMEDIR = os.path.expandvars("$HOME")
        # self.WORKDIR = os.path.expandvars("$WORK")
        self.HOMEDIR = '/n/conroyfs1/pac/FAL'
        self.WORKDIR = '/n/conroyfs1/pac/FAL'
        self.exedir = self.HOMEDIR+"/bin/"
        self.datadir = self.HOMEDIR+"/data/"
        self.bigdatadir = self.WORKDIR+"/data/"

        # create the glue  
        self.glue = FALGlue.glue()

        # Call selmod to set model atm and other star specific parameters
        (self.synbegvar,self.instparstr,self.modatm) = selmod(self.starpars)

        # move mod atm into memory
        self._fastcopy(self.modatm,'/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
        self.atmomod = '/dev/shm/FAL/{0}/mod.dat'.format(self.ID)

        # set up some useful strings
        self.rotatevar = ("    1\n{VROT:4.3f}")
        self.macpar = "MACRO     {MACVEL:3.1f}       KM                  COMMENT FIELD"

        # move some static files into memory
        # MOLECULES
        try:
            os.remove('/dev/shm/FAL/{0}/molecules.dat'.format(self.ID))
        except OSError:
            pass
        # BOB's MOLECULES FILE
        self._fastcopy(self.datadir+'molecules.dat','/dev/shm/FAL/{0}/molecules.dat'.format(self.ID))
        self.molecules = '/dev/shm/FAL/{0}/molecules.dat'.format(self.ID)
        
        # CONTINUA
        try:
            os.remove('/dev/shm/FAL/{0}/continua.dat'.format(self.ID))
        except OSError:
            pass
        # BOB's CONTINUA FILE
        self._fastcopy(self.datadir+'continua.dat','/dev/shm/FAL/{0}/continua.dat'.format(self.ID))
        self.continua = '/dev/shm/FAL/{0}/continua.dat'.format(self.ID)
        
        # He1Tables
        try:
            os.remove('/dev/shm/FAL/{0}/he1tables.dat'.format(self.ID))
        except OSError:
            pass
        self._fastcopy(self.datadir+'he1tables.dat','/dev/shm/FAL/{0}/he1tables.dat'.format(self.ID))
        self.he1tables = '/dev/shm/FAL/{0}/he1tables.dat'.format(self.ID)

        # SPECTRV.INPUT
        try:
            os.remove('/dev/shm/FAL/{0}/spectrv.input'.format(self.ID))
        except OSError:
            pass
        self._fastcopy(self.datadir+'spectrv.input','/dev/shm/FAL/{0}/spectrv.input'.format(self.ID))
        self.spectrvin = '/dev/shm/FAL/{0}/spectrv.input'.format(self.ID)

    def __delete_(self):
        pass

    def _callpro(self,function,inputstr=None,inpipe=None,verbose=None):
        """
        general function to call fortran code
        """
        _inputstr = inputstr
        # set verbose
        if verbose == None:
            _FNULL = self.FNULL
        elif verbose == True:
            _FNULL = sys.stdout
        elif verbose == False:
            _FNULL = open(os.devnull, 'w')
        else:
            _FNULL = self.FNULL

        if inpipe != None:
            pro = subprocess.Popen([self.exedir+function+".exe","_"+self.ID],
                                   stdin=open(inpipe,'r'),stdout=_FNULL)
        else:
            pro = subprocess.Popen([self.exedir+function+".exe","_"+self.ID],
                                   stdin=subprocess.PIPE,stdout=_FNULL)

        output = pro.communicate(_inputstr)
        return output

    def _stringIOafy(self,strin):
        """
        general function to input in string into StringIO object
        """
        output = cStringIO.StringIO()
        output.write(strin)
        output.close()
        return output

    def _fastcopy(self,src,dst,buffer_size=-1):#10485760):
        """
        Function that does a fast copy using buffers and binary fmt
        """
        # Optimize buffer for small files
        # buffer_size = min(buffer_size,os.path.getsize(src))
        # if(buffer_size == 0):
        #     buffer_size = 1024
        try:
            with open(src,'rb') as fsrc:
                with open(dst,'wb') as fdst:
                    shutil.copyfileobj(fsrc,fdst,buffer_size)
        except TypeError:
            shutil.copy(src,dst)
            
    def _makesym(self,src,outname):
        """
        Function that takes a file, copies to memory, and sets up a symlink
        """
        self._fastcopy(src,'/dev/shm/FAL/{0}/{1}'.format(self.ID,outname))
        os.symlink('/dev/shm/FAL/{0}/{1}'.format(self.ID,outname),'./'+outname)

    def _mvsym(self,src,outname):
        """
        Function that takes a file, copies to memory, and sets up a symlink
        """
        self._fastcopy('/dev/shm/FAL/{0}/{1}'.format(self.ID,src),
                       '/dev/shm/FAL/{0}/{1}'.format(self.ID,outname))
        os.remove('/dev/shm/FAL/{0}/{1}'.format(self.ID,src))
        os.unlink(src)
        os.symlink('/dev/shm/FAL/{0}/{1}'.format(self.ID,outname),'./'+outname)

    def _rmsym(self,src,verbose=True):
        """
        Function that removes sym and file in memory
        """
        try:
            os.remove('/dev/shm/FAL/{0}/{1}'.format(self.ID,src))
        except OSError:
            if verbose:
                print('WARNING: Could not find /dev/shm/FAL/{0}/{1}'.format(self.ID,src))
        try:
            os.unlink(src)
        except OSError:
            if verbose:
                print('WARNING: Could not find {0}'.format(src))
        
    def _fileprep(self,filesdict):
        """
        function to make sure all infiles exisit and new outfiles don't
        """

        # create filename strings
        localfilenamestr = "./{NAME}"
        memfilenamestr = "/dev/shm/FAL/{ID}/{NAME}"
        
        # check exisitance of infiles
        if len(filesdict["infiles"].keys()) > 0:
            for fname in filesdict['infiles'].keys():
                # create filename strings
                lfs = localfilenamestr.format(NAME=fname)
                mfs = memfilenamestr.format(NAME=fname,ID=self.ID)
                
                # check local file
                if not os.path.isfile(lfs):
                    print("Error: Missing Input file -> {NAME}".format(NAME=lfs))
                    return "Error: Missing Input file -> {NAME}".format(NAME=lfs)

                # check file in memory
                if not os.path.isfile(mfs):
                    print("Error: Missing Input file -> {NAME}".format(NAME=mfs))
                    return "Error: Missing Input file -> {NAME}".format(NAME=mfs)
        
        # check exisitance of write-to files
        if len(filesdict["writeto"].keys()) > 0:
            for fname in filesdict['writeto'].keys():
                # create filename strings
                lfs = localfilenamestr.format(NAME=fname)
                mfs = memfilenamestr.format(NAME=fname,ID=self.ID)
                
                # check local file
                if not os.path.isfile(lfs):
                    return "Error: Missing Write-to file -> {NAME}".format(NAME=lfs)

                # check file in memory
                if not os.path.isfile(mfs):
                    return "Error: Missing Write-to file -> {NAME}".format(NAME=mfs)

        # create new symbolc links and files for output
        if len(filesdict["newfiles"].keys()) > 0:
            for fname in filesdict['newfiles'].keys():
                # create filename strings
                lfs = localfilenamestr.format(NAME=fname)
                mfs = memfilenamestr.format(NAME=fname,ID=self.ID)

                # check file in memory
                try:
                    os.remove(mfs)
                except OSError:
                    pass
                if filesdict['newfiles'][fname] == 'bin':
                    #ff = spIO.FortranFile(mfs,'w')
                    ff = open(mfs,'wb')
                    f = np.fromfile(ff,count=0)
                    
                elif filesdict['newfiles'][fname] == 'ascii':
                    ff = open(mfs,'w')
                else:
                    return "Error: Did not understand file type on -> {NAME}".format(NAME=mfs)
                ff.close()

                # check local file
                try:
                    os.unlink(lfs)
                except OSError:
                    pass
                os.symlink(mfs,lfs)
        
    def cleandir(self):
        """
        Function to clean all fort* symlinks and /dev/shm/FAL/ID/fort* files
        """
        # first do symlinks
        symlinklist = glob.glob('./fort*')
        for ff in symlinklist:
            os.unlink(ff)
        
        # then do files in memory
        filelist = glob.glob('/dev/shm/FAL/{0}/fort*'.format(self.ID))
        for ff in filelist:
            os.remove(ff)

    def synbeg(self,indict,clobber=False,verbose=None):
        """
        Run SYNBEG code

        Reads In: 
            command line

        Writes Into:
            None

        New Out: 
            fort.10 (bin)
            fort.12 (bin)
            fort.14 (bin)
            fort.19 (bin)
            fort.20 (bin)
            fort.93 (bin)
        """
        # clobber old run if wanted
        if clobber == True:
            # first delete any *_* files in memory
            memfiles = glob.glob('/dev/shm/FAL/{0}/fort.*'.format(self.ID))
            for mf in memfiles:
                os.unlink(mf)
            # clear local directory of broken links
            locfiles = glob.glob('./*')
            for lf in locfiles:
                try:
                    os.stat(lf)
                except OSError:
                    os.unlink(lf)

        # read in information from input dictionary
        self.WSTART = self.starpars['WSTART']
        self.WEND = self.starpars['WEND']
        if 'RESOL' in self.starpars.keys():
            self.RESOL = self.starpars['RESOL']
        else:
            self.RESOL = 3000000.0
        if 'PRED' in self.starpars.keys():
            self.PRED = self.starpars['PRED']
        else:
            self.PRED = 1
        if 'LINOUT' in self.starpars.keys():
            self.LINOUT = self.starpars['LINOUT']
        else:
            self.LINOUT = 30
        if 'TOL' in self.starpars.keys():
            self.TOL = self.starpars['TOL']
        else:
            self.TOL = 1e-3
        if 'OUTRES' in self.starpars.keys():
            self.OUTRES = self.starpars['OUTRES']
        else:
            self.OUTRES = None
        
        # remove any extant RUNINFO.dat
        try:
            os.remove('./RUNINFO.dat')
        except OSError:
            pass

        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {}
        filesdict['writeto'] = {}
        filesdict['newfiles'] = ({'fort.10':'bin','fort.12':'bin',
                                  'fort.14':'bin','fort.19':'bin',
                                  'fort.20':'bin','fort.93':'bin'})
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")
            

        # write in information into input string
        if verbose:
            print("Running SynBeg")
        synbegvar_i = self.synbegvar.format(WSTART=self.WSTART,WEND=self.WEND,RESOL=self.RESOL,PRED=self.PRED,TOL=self.TOL,LINOUT=self.LINOUT,OUTRES=self.OUTRES)
        self.synbegout = self._callpro("synbeg",synbegvar_i,verbose=verbose)
        if verbose:
            print("Finished SynBeg")

        return (self.synbegout,self.ID)

    def readlines(self,rtype=None,rlinedict=None,verbose=None,MASTERLL=None):
        """
        function that determines if you should use a punch500 or 
        read in all the individual lines
        """

        if (rtype == "readlast"):
            self.rgfalldel(verbose=verbose,reuse_ll=True)
        elif (rtype == 'readmaster'):
            self.rmaster(MASTERLL=MASTERLL,verbose=verbose)
        elif (rtype == 'readall'):
            self.rlinefunc(rlinedict,verbose=verbose)
        else:
            self.rgfalldel(verbose=verbose,reuse_ll=False,userll=rtype)

    def rmaster(self,MASTERLL=None,verbose=None):
        """
        Run RPUNCHBIN code

        Reads In:
            fort.11 (bin) [master line list]

        Writes Into:
            fort.12 (bin)
            fort.14 (bin)
            fort.19 (bin)
            fort.20 (bin)
            fort.93 (bin)

        New Out: 
            None
        """
        # Master lines file is gigantic, so don't copy into memory just sym link it
        if MASTERLL == None:
            MASTERLL = (['/n/conroyfs1/pac/FAL/data/LL/KuruczLL_1400_1900.bin',
                '/n/conroyfs1/pac/FAL/data/LL/KuruczH2OLL_1400_1900.bin'])

        # for each line list in MASTERLL, run rpunchbin
        for MLL in MASTERLL:

            # link master line list file
            if os.path.isfile("fort.11"):
                self._rmsym('fort.11',verbose=verbose)
            os.symlink(MLL,'fort.11')

            # run read master list program
            self.rmasterout = self._callpro("rpunchbin",verbose=verbose)


        # move fort.14 to fort.99 and clean up
        self._mvsym('fort.14','fort.99')
        
        return (self.rmasterout,self.ID)


    def rgfalldel(self,verbose=None,reuse_ll=True,userll=None):
        """
        Run RMOLECASC codes
        
        Reads In: 
            fort.11 (ascii)[Line List]

        Writes Into:
            fort.12 (bin)
            fort.14 (bin)
            fort.19 (bin)
            fort.20 (bin)
            fort.93 (bin)

        New Out: 
            None
        """
        if reuse_ll:
            pass
        else:
            # write line list into fort.11
            if os.path.isfile("fort.11"):
                self._rmsym('fort.11',verbose=verbose)
            # check to see if user line list is np table ll or an ascii path
            if type(userll).__name__ == 'Table':
                self.glue.writenp(userll,'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
                os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')
            else:
                self._makesym(userll,'fort.11')
        
        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.11':'ascii'}
        filesdict['writeto'] = ({'fort.12':'bin','fort.14':'bin','fort.19':'bin',
                                  'fort.20':'bin','fort.93':'bin'})
        filesdict['newfiles'] = {}
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        #print("Running RgfAlldel")
        # self.rgfalldelout = self._callpro("rgfalldel",verbose=verbose)
        self.rgfalldelout = self._callpro("rpunchxx",verbose=verbose)
        #print("Finished RgfAlldel")

        # move fort.14 to fort.99 and clean up
        self._mvsym('fort.14','fort.99')
        
        return (self.rgfalldelout,self.ID)


    def rlinefunc(self,rlinedict={"atoms":True,"moles":True},verbose=None):
        """
        Run read line code to manually read in all data
        """

        try:
            if rlinedict['atoms']==True:
                # do atomic lines
                self.ratomic(verbose=verbose)
        except KeyError:
            pass
        try:
            if rlinedict['predict']==True:
                # do predicted lines
                self.readpredlines(verbose=verbose)
        except KeyError:
            pass
        try:
            if rlinedict['moles']==True:
                self.molefiles = ([
                    'alopatrascu.asc', # AlO
                    'nah.dat', # NaH
                    'voax.asc','vobx.asc','vocx.asc', # VO
                    'fehfx.dat', # FeH
                    'h2bx.dat','h2cx.dat','h2xx.dat', # H2
                    'hdxx.dat', #HD
                    'mgh.dat', #MgH
                    # 'mghax.dat','mghbx.dat', #MgH
                    'mgodaily.dat', #MgO
                    'nhax.dat','nhca.dat', # NH
                    'chjorg.dat', #CH
                    # 'cnaxbrooke.dat','cnbxbrooke.dat','cnxx12brooke.dat', # CN
                    'cnax.dat','cnbx.dat', #CN
                    'c2ax.dat','c2ba.dat','c2da.dat','c2ea.dat', #C2
                    'coax.dat','coxx.dat', # CO
                    # 'ohax.dat','ohxx.dat', # OH
                    'ohaxupdate.asc','ohxxgoldman.asc',# OH
                    'sihax.dat', # SiH
                    'sioax.dat','sioex.dat','sioxx.dat', #SiO
                    'crhax.dat', # CrH
                    'cah.dat', # CaH
                    ])
                # read molecular files
                self.rmolecascout = {}
                for mf in self.molefiles:
                    try:
                        self.readmol(mf,verbose=verbose)
                    except IOError:
                        print('!!!!! Missing {0} !!!!!'.format(mf))
                        raise
        except KeyError:
            pass
        try:
            if rlinedict['TiO']==True:
                # do TiO lines
                self.readmol_TiO(verbose=verbose)
            if os.path.isfile("fort.48"):
                os.unlink('fort.48')
            if os.path.isfile('/dev/shm/FAL/{0}/fort.48'.format(self.ID)):
                os.remove('/dev/shm/FAL/{0}/fort.48'.format(self.ID))
            if os.path.isfile("fort.11"):
                os.unlink('fort.11')
            if os.path.isfile('/dev/shm/FAL/{0}/fort.11'.format(self.ID)):
                os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))
        except KeyError:
            pass
        try:
            if rlinedict['H2O']==True:
                # do H20 lines
                self.readmol_H2O(verbose=verbose)
            if os.path.isfile("fort.48"):
                os.unlink('fort.48')
            if os.path.isfile('/dev/shm/FAL/{0}/fort.48'.format(self.ID)):
                os.remove('/dev/shm/FAL/{0}/fort.48'.format(self.ID))
            if os.path.isfile("fort.11"):
                os.unlink('fort.11')
            if os.path.isfile('/dev/shm/FAL/{0}/fort.11'.format(self.ID)):
                os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))
        except KeyError:
            pass
        
        try:
            if rlinedict['injectlines'] != None:
                # inject lines into line list
                self.rinjectlines(verbose=verbose,injectll=rlinedict['injectlines'])
        except KeyError:
            pass

        # move fort.14 to fort.99
        self._mvsym('fort.14','fort.99')

        return ("COMPLETED ALL CALLS IN RLINE",self.ID)

    def readmol_TiO(self,verbose=None):
        """
        helper function to do molecule TiO file read in
        """

        if verbose:
            print("Running RSCHWENK")

        # write TiO line list files into fort.11 and fort.48
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        os.symlink(self.bigdatadir+'MOLECULES/schwenke.bin','fort.11')
        os.symlink(self.bigdatadir+'MOLECULES/schwenke.bin','/dev/shm/FAL/{0}/fort.11'.format(self.ID))

        if os.path.isfile("fort.48"):
            self._rmsym('fort.48',verbose=verbose)
        os.symlink(self.bigdatadir+'MOLECULES/eschwenke.bin','fort.48')
        os.symlink(self.bigdatadir+'MOLECULES/eschwenke.bin','/dev/shm/FAL/{0}/fort.48'.format(self.ID))        

        self.rmoleout_tio = self._callpro("rschwenk",verbose=verbose)

        if os.path.isfile("fort.48"):
            os.unlink('fort.48')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.48'.format(self.ID)):
            os.remove('/dev/shm/FAL/{0}/fort.48'.format(self.ID))
        if verbose:
            print("Finished RSCHWENK")
    
    def readmol_H2O(self,verbose=None):
        """
        helper function to do molecule H2O file read in
        """

        if verbose:
            print("Running RH2OSLOW")

        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        os.symlink(self.bigdatadir+'MOLECULES/h2oslowfix.bin','fort.11')
        os.symlink(self.bigdatadir+'MOLECULES/h2oslowfix.bin','/dev/shm/FAL/{0}/fort.11'.format(self.ID))        

        if os.path.isfile("fort.48"):
            self._rmsym('fort.48',verbose=verbose)
        os.symlink(self.bigdatadir+'MOLECULES/eh2opartridge.bin','fort.48')
        os.symlink(self.bigdatadir+'MOLECULES/eh2opartridge.bin','/dev/shm/FAL/{0}/fort.48'.format(self.ID))        

        self.rmoleout_h2o = self._callpro("rh2oslow",verbose=verbose)

        if os.path.isfile("fort.48"):
            os.unlink('fort.48')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.48'.format(self.ID)):
            os.remove('/dev/shm/FAL/{0}/fort.48'.format(self.ID))

        if verbose:
            print("Finished RH2OSLOW")

    def readmol(self,molfile,verbose=None):
        """
        helper function to do molecule file read in
        """
        # write molfile line list into fort.11
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        os.symlink(self.bigdatadir+'MOLECULES/{0}'.format(molfile),'fort.11')
        os.symlink(self.bigdatadir+'MOLECULES/{0}'.format(molfile),'/dev/shm/FAL/{0}/fort.11'.format(self.ID))        
        
        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.11':'ascii'}
        filesdict['writeto'] = ({'fort.12':'bin','fort.14':'bin','fort.19':'bin',
            'fort.20':'bin','fort.93':'bin'})
        filesdict['newfiles'] = {}
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        if verbose:
            print("Running RMolecASC on {0}".format(molfile))
        self.rmolecascout[molfile] = self._callpro("rmolecasc",verbose=verbose)
        if verbose:
            print("Finished RMolecASC on {0}".format(molfile))

    def readpredlines(self,verbose=None):
        """
        helper function to read in predicted lines
        """
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)

        # Predicted lines file is gigantic, so don't copy into memory just sym link it
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        # os.symlink(self.bigdatadir+'/gfpred27sep15.bin','fort.11')
        os.symlink(self.bigdatadir+'/gfpred29dec2014.bin','fort.11')
        self.rmoleout = self._callpro("rpredict",verbose=verbose)

    def ratomic(self,verbose=None):
        """
        helper function to read in atomic lines
        """
        # write atomic line list into fort.11
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        # os.symlink(self.bigdatadir+'/gfall18feb16.dat','fort.11')
        # os.symlink(self.bigdatadir+'/gfall18feb16.dat','/dev/shm/FAL/{0}/fort.11'.format(self.ID))
        # os.symlink(self.bigdatadir+'/gfall05jun16.dat','fort.11')
        # os.symlink(self.bigdatadir+'/gfall05jun16.dat','/dev/shm/FAL/{0}/fort.11'.format(self.ID))
        os.symlink(self.bigdatadir+'/gfall18feb16.dat','fort.11')
        os.symlink(self.bigdatadir+'/gfall18feb16.dat','/dev/shm/FAL/{0}/fort.11'.format(self.ID))


        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.11':'ascii'}
        filesdict['writeto'] = ({'fort.12':'bin','fort.14':'bin','fort.19':'bin',
            'fort.20':'bin','fort.93':'bin'})
        filesdict['newfiles'] = {}
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        # print("Running RGFALL")
        self.ratomicout = self._callpro("rgfall",verbose=verbose)
        # print("Finished RGFALL")

    def rinjectlines(self,verbose=None,injectll=None):
        """
        Run RPUNCHXX ON INJECTED LINES
        
        Reads In: 
            fort.11 (ascii)[Line List]

        Writes Into:
            fort.12 (bin)
            fort.14 (bin)
            fort.19 (bin)
            fort.20 (bin)
            fort.93 (bin)

        New Out: 
            None
        """
        # write line list into fort.11
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        self._makesym(injectll,'fort.11')
        
        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.11':'ascii'}
        filesdict['writeto'] = ({'fort.12':'bin','fort.14':'bin','fort.19':'bin',
                                  'fort.20':'bin','fort.93':'bin'})
        filesdict['newfiles'] = {}
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        #print("Running RgfAlldel")
        # self.rgfalldelout = self._callpro("rgfalldel",verbose=verbose)
        self.rinjectllout = self._callpro("rpunchxx",verbose=verbose)
        #print("Finished RgfAlldel")


    def xnfpelsyn(self,verbose=None):
        """
        Run XNFPELSYN code
        
        Reads In: 
            fort.2  (ascii)[molecules]
            fort.17 (bin)[continua]
            fort.18 (ascii)[he1lines]

        Writes Into:
            fort.10 (bin)
            
        New Out: 
            fort.35 (ascii)
        """
        # write links to molecules, continua, and he1tables
        self._makesym(self.molecules,'fort.2')
        self._makesym(self.continua,'fort.17')
        self._makesym(self.he1tables,'fort.18')
        
        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.2':'ascii','fort.17':'ascii','fort.18':'ascii'}
        filesdict['writeto'] = ({'fort.10':'bin'})
        filesdict['newfiles'] = {'fort.35':'ascii'}
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        # print("Running xnfpelsyn")
        self.xnfpelsynout = self._callpro("xnfpelsyn",
                                          inpipe="{0}".format(self.atmomod),
                                          verbose=verbose)
        # print("Finished xnfpelsyn")
        
        return (self.xnfpelsynout,self.ID)

    def syn(self,verbose=None,speed=None):
        """
        Run SYNTHE code

        Reads In: 
           fort.93
           fort.10
        Writes Into:
           fort.12
           fort.19
           fort.20
           fort.99
        New Out: 
           fort.7
           fort.8
           fort.9
           fort.13
           fort.14
           fort.15
        Delete At Completion
           fort.7
           fort.13
           fort.14
           fort.15

        """

        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.93':'bin','fort.10':'bin'}
        filesdict['writeto'] = {'fort.12':'bin','fort.19':'bin','fort.20':'bin','fort.99':'bin'}
        filesdict['newfiles'] = ({'fort.7':'bin','fort.8':'bin',
                                  'fort.9':'bin','fort.13':'bin',
                                  'fort.14':'bin','fort.15':'bin'})
        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")
        
        # write in information into input string
        # print("Running Synthe")
        if speed == 'slow':
            self.synout = self._callpro("synthe_slow",verbose=verbose)
        elif speed == 'fast':
            self.synout = self._callpro("synthe_fast",verbose=verbose)
        elif speed == 'ultrafast':
            self.synout = self._callpro("synthe_ultrafast",verbose=verbose)
        else:
            self.synout = self._callpro("synthe",verbose=verbose)
        # print("Finished Synthe")

        self._rmsym('fort.7',verbose=False)
        self._rmsym('fort.13',verbose=False)
        self._rmsym('fort.14',verbose=False)
        self._rmsym('fort.15',verbose=False)
                
        return (self.synout,self.ID)

    def spectrv(self,verbose=None,tau=False):
        """
        Run SPECTRV code

        Reads In: 
            fort.5 [model atm]
            fort.25 [spectrv.input]
            fort.9
            fort.10
        Writes Into:
        
        New Out: 
            fort.7 (bin)
            fort.16 (ascii)
            fort.33 (ascii)
            fort.20 (bin)
        Delete At Completion
            fort.20
            fort.9
        """
        # make mod atm link to fort.5 and spectrv.input as fort.25
        self._makesym(self.atmomod,'fort.5')
        self._makesym(self.spectrvin,'fort.25')
        
        if os.path.exists('./LINEINFO.dat'):
            os.remove('./LINEINFO.dat')
     
        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.5':'ascii','fort.25':'ascii','fort.9':'bin','fort.25':'bin'}
        filesdict['writeto'] = {}
        filesdict['newfiles'] = ({'fort.7':'bin','fort.16':'ascii',
                                  'fort.33':'ascii','fort.20':'bin'})

        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")
            

        # write in information into input string
        # print("Running SpectrV")
        if tau:
            self.spectrvout = self._callpro("spectrv_tau",verbose=verbose)
        else:
            self.spectrvout = self._callpro("spectrv",verbose=verbose)
        # print("Finished SpectrV")
        
        if os.path.isfile('fort.9'):
            os.unlink('fort.9')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.9'.format(self.ID)):
            os.remove('/dev/shm/FAL/{0}/fort.9'.format(self.ID))
        if os.path.isfile('fort.20'):
            os.unlink('fort.20')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.20'.format(self.ID)):
            os.remove('/dev/shm/FAL/{0}/fort.20'.format(self.ID))
        
        self._mvsym('fort.7','fort.sol')
        if tau:
            shutil.copy('/dev/shm/FAL/{0}/fort.33'.format(self.ID),'ROSS.tau')

        return (self.spectrvout,self.ID)

    def rotate(self,VROT=0.0,verbose=None):
        """
        Run ROTATE code

        Reads In: 
             fort.1
        Writes Into:
             fort.19
        New Out: 
             ROTX (bin) X = # rotation velocities
        Delete At Completion:
        """
        # link fort.sol to fort.1
        self._makesym('/dev/shm/FAL/{0}/fort.sol'.format(self.ID),'fort.1')
                
        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.1':'bin'}
        filesdict['writeto'] = {'fort.19':'bin'}
        filesdict['newfiles'] = {'ROT1':'bin'}

        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")
            

        # write in information into input string
        # print("Running Rotate")
        rotatestr = self.rotatevar.format(VROT=VROT)
        self.rotateout = self._callpro("rotate",rotatestr,verbose=verbose)
        # print("Finished Rotate")

        return (self.rotateout,self.ID)
    
    def broaden(self,inspec,VMAC=0.0,broadtype=None,WLreg=None,write=False,verbose=None):
        if broadtype==None:
            raise ValueError('Must define broadening type')

        # clean out fort.21 (input)
        if os.path.exists('fort.21'):
            os.unlink('fort.21')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.21'):
            os.remove('/dev/shm/FAL/{0}/fort.21')
                    
        # mv input file into fort.21
        self._makesym('/dev/shm/FAL/{0}/{1}'.format(self.ID,inspec),'fort.21')

        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.21':'bin'}
        filesdict['writeto'] = {}
        filesdict['newfiles'] = {'fort.22':'bin'}

        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        # now decide which type of broadening and run
        if broadtype=="MAC":
            tag = 'mac'
            macstr = self.macpar.format(MACVEL=VMAC)
            print(macstr)
            self.broadout = self._callpro(
                "broadenx",inputstr=macstr,
                verbose=verbose)

        elif broadtype=="INSTRUMENT":
            tag = 'inst'
            # determine which instrument settings should be used
            if WLreg=='OPT':
                intpars = self.instparstr['OPT']
            elif WLreg=='HBAND':
                intpars = self.instparstr['HBAND']
            else:
                intpars = self.instparstr['HBAND']

            parset = intpars.keys()
            if len(intpars.keys()) == 1:
                # the case with only one instrumental broadening (likely just a gaussian)
                self.broadout = self._callpro("broadenx",inputstr=intpars[parset[0]],verbose=verbose)
            else:
                # the special case of the solar line profile with a SINC and a gaussian
                intpars_1 = intpars[parset[0]]
                intpars_2 = intpars[parset[1]]

                # first run sinc
                self.broadout1 = self._callpro("broadenx",inputstr=intpars_1,verbose=verbose)
                # mv fort.22 into fort.21 as into into gauss
                os.unlink('fort.21')
                os.remove('/dev/shm/FAL/{0}/fort.21'.format(self.ID))
                self._mvsym('fort.22','fort.21')
                ff = open('/dev/shm/FAL/{0}/fort.22'.format(self.ID),'wb')
                f = np.fromfile(ff,count=0)
                ff.close()
                os.symlink('/dev/shm/FAL/{0}/fort.22'.format(self.ID),'fort.22')
                self.broadout2 = self._callpro("broadenx",inputstr=intpars_2,verbose=verbose)
        else:
            raise ValueError("Did not understand broadening type (MAC/INSTRUMENT)")
        
        # move final fort.22 into an output file with the right tag
        if os.path.exists(inspec+'_{0}'.format(tag)):
            os.remove(inspec+'_{0}'.format(tag))
        self._mvsym('fort.22',inspec+'_{0}'.format(tag))

    def transmit(self):
        pass

    def writespec(self,inspec,verbose=None):
        """
        Function to write out binary spectrum to ascii
        """
        
        # remove any old ascii output files
        outfiles = ['headinfo.dat','lineinfo.dat','specfile.dat']
        for ofile in outfiles:
            try:
                os.remove(ofile)
            except OSError:
                pass

        # remove any old fort.1 and copy binary file to fort.1
        if os.path.exists('fort.1'):
            os.unlink('fort.1')
        shutil.copy(inspec,'fort.1')

        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.1':'bin'}
        filesdict['writeto'] = {}
        filesdict['newfiles'] = {}

        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        self.writespecout = self._callpro('syntoascanga',verbose=verbose)

        
    def at12tosyn(self):
        pass

    def archive(self):
        """
        Archive Button
        """
        # make copy of files in INT directory so that everything can be restarted fresh
        if not os.path.exists('/dev/shm/FAL/{0}/INT'.format(self.ID)):
            os.makedirs('/dev/shm/FAL/{0}/INT'.format(self.ID))
        fortlist = ['fort.10','fort.12','fort.14','fort.19','fort.20','fort.93']
        for ft in fortlist:
            self._fastcopy('/dev/shm/FAL/{0}/{1}'.format(self.ID,ft),
                           '/dev/shm/FAL/{0}/INT/{1}'.format(self.ID,ft))

    def reset(self):
        """
        RESET BUTTON
        """
        # remove all all old files except fort.11
        fortlist = glob.glob('fort.*') + glob.glob('ROT*')
        for ft in fortlist:
            self._rmsym(ft)
        if os.path.exists('LINEINFO.dat'):
            os.remove('LINEINFO.dat')

        # copy files for INT back into memory and link them
        fortlist = ['fort.11','fort.10','fort.12','fort.14','fort.19','fort.20','fort.93']
        for ft in fortlist:
            self._fastcopy('/dev/shm/FAL/{0}/INT/{1}'.format(self.ID,ft),
                '/dev/shm/FAL/{0}/{1}'.format(self.ID,ft))
            os.symlink('/dev/shm/FAL/{0}/{1}'.format(self.ID,ft),ft)

    @property
    def atmomod(self):
        """
        This is a cached version of `get_atmomod`, which you can
        override to set the atmomod. If that method is not overridden, 
        this method always returns None.

        We lazy-load the atmomod the first time this method is called
        and cache the result after that.
        """
        if not hasattr(self, "_atmomod"):
            self._atmomod = self.get_atmomod()
        return self._atmomod

    @atmomod.setter
    def atmomod(self,value):
        self._atmomod = value

    def get_atmomod(self):
        """Override to determine the atmomod"""
        return None

    @property
    def linelist(self):
        """
        This is a cached version of `get_linelist`, which you can
        override to set the linelist. If that method is not overridden, 
        this method always returns None.

        We lazy-load the linelist the first time this method is called
        and cache the result after that.
        """
        if not hasattr(self, "_linelist"):
            self._linelist = self.get_linelist()
        return self._linelist

    @linelist.setter
    def linelist(self,value):
        self._linelist = value

    def get_linelist(self):
        """Override to determine the linelist"""
        return None

if __name__ == "__main__":
    indict = {'WSTART':516.4,'WEND':519.0,'RESOL':2000000.0,'VROT':0.0,'MACVEL':1.5,'PRED':0}
    syn = synthe(verbose=False,clobber=True)
    syn(indict)
