
from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import os,sys,glob,time,filecmp
import subprocess
import cStringIO
import shutil
from scipy import io as spIO
import numpy as np

__all__ = ["synthepkg"]


class synthe(object):

    def __init__(self,ID=None,verbose=False,clobber=False,rline=False,LL=False,rlinedict={},arct_bool=False):
        super(synthe,self).__init__()
        
        # set if working with Arcturus
        self.arct_bool = arct_bool

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
            memfiles = glob.glob('/dev/shm/FAL/{0}/*_*'.format(self.ID))
            for mf in memfiles:
                os.unlink(mf)
            # clear local directory of broken links
            locfiles = glob.glob('./*')
            for lf in locfiles:
                try:
                    os.stat(lf)
                except OSError:
                    os.unlink(lf)
            
        # Set up exec directory and data directory
        HOMEDIR = os.path.expandvars("$HOME")
        self.exedir = HOMEDIR+"/FAL/PYTHON/bin/"
        self.datadir = HOMEDIR+"/FAL/PYTHON/data/"
        self.moldatadir = '/work/02349/cargilpa/FAL/DATA/MOLECULES/'

        # define if writes to /dev/null or StdOUT
        if verbose == True:
            self.FNULL = sys.stdout
        elif verbose == False:
            self.FNULL = open(os.devnull, 'w')
        else:
            raise ValueError('Do not understand verbose flag')
            
        # copy input line list and model atm into memory
        # then set up symbolic links

        # LINELIST ->

        # RLINE MEANS TO RUN A NEW LINE LIST, NOT PUNCH500
        if rline:
            self.molefiles = ([
                'alopatrascu.asc', # AlO
                'nah.dat', # NaH
                'voax.asc','vobx.asc','vocx.asc', # VO
                'fehfx.dat', # FeH
                'h2bx.dat','h2cx.dat','h2xx.dat', # H2
                'hdxx.dat', #HD
                'mgh.dat', #MgH
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

            self.rline = True
            self.savetomem = False
            if self.savetomem:
                if "atoms" in rlinedict.keys():
                    # RLINE ATOMIC LINE LIST FILE
                    src1ll = '/work/02349/cargilpa/FAL/DATA/gfall25jul15.dat'
                    src2ll = '/dev/shm/FAL/{0}/gfall25jul15.dat'.format(self.ID)
                    distll = '/dev/shm/FAL/{0}/atomicll.dat'.format(self.ID)

                    if os.path.exists(src2ll):
                        if filecmp.cmp(src1ll,src2ll):
                            # print('... using old {0}'.format(src2ll))
                            pass
                        else:
                            # print('... copying {0}'.format(src1ll))
                            self._fastcopy(src1ll,src2ll)
                    else:
                        # print('... copying {0}'.format(src1ll))
                        self._fastcopy(src1ll,src2ll)

                    self._fastcopy(src2ll,distll)
                    self.atomiclinelist = distll

                if "moles" in rlinedict.keys():
                    # COPY THE MOLECULE FILES INTO MEMORY
                    for molf in self.molefiles:
                        srcpath = self.moldatadir+'{0}'.format(molf)
                        dstpath = '/dev/shm/FAL/{0}/{1}'.format(self.ID,molf)
                    
                        if os.path.isfile(dstpath):
                            if filecmp.cmp(srcpath,dstpath):
                                # print('... using old {0}'.format(dstpath))
                                pass
                            else:   
                                # print('... copying {0}'.format(srcpath))
                                self._fastcopy(srcpath,dstpath)
                        else:
                            # print('... copying {0}'.format(srcpath))
                            self._fastcopy(srcpath,dstpath)
                if "TiO" in rlinedict.keys():
                    # COPY THE TiO FILES INTO MEMORY
                    molefiles = (['schwenke.bin','eschwenke.bin'])
                    for molf in molefiles:
                        srcpath = self.moldatadir+'{0}'.format(molf)
                        dstpath = '/dev/shm/FAL/{0}/{1}'.format(self.ID,molf)
                        if os.path.isfile(dstpath):
                            if filecmp.cmp(srcpath,dstpath):
                                # print('... using old {0}'.format(dstpath))
                                pass
                            else:
                                # print('... copying {0}'.format(srcpath))
                                self._fastcopy(srcpath,dstpath)
                        else:
                            # print('... copying {0}'.format(srcpath))
                            self._fastcopy(srcpath,dstpath)
                if "H2O" in rlinedict.keys():
                    # COPY THE H2O FILES INTO MEMORY
                    molefiles = (['h2oslowfix.bin','eh2opartridge.bin'])
                    for molf in molefiles:
                        srcpath = self.moldatadir+'{0}'.format(molf)
                        dstpath = '/dev/shm/FAL/{0}/{1}'.format(self.ID,molf)
                        if os.path.isfile(dstpath):
                            if filecmp.cmp(srcpath,dstpath):
                                # print('... using old {0}'.format(dstpath))
                                pass
                            else:
                                # print('... copying {0}'.format(srcpath))
                                self._fastcopy(srcpath,dstpath)
                        else:
                            # print('... copying {0}'.format(srcpath))
                            self._fastcopy(srcpath,dstpath)
        else:
            if LL:
                self.linelist = '/dev/shm/FAL/{0}/ll.dat'.format(self.ID)
            else:
                try:
                    os.remove('/dev/shm/FAL/{0}/ll.dat'.format(self.ID))
                except OSError:
                    pass
                if rlinedict["ll"] == 'punch500':
                    # BOB's FULL PUNCH500 FILE
                    self._fastcopy(self.datadir+'punch500_new.del','/dev/shm/FAL/{0}/ll.dat'.format(self.ID))
                    # Bob's PUNCH500A FILE
                    # self._fastcopy(self.datadir+'punch500a.del','/dev/shm/FAL/ll.dat_{0}'.format(self.ID)) 
                    # BOB's PUNCH500B FILE
                    # self._fastcopy(self.datadir+'punch500b.del','/dev/shm/FAL/ll.dat_{0}'.format(self.ID))
                    # BOB's PUNCH500C FILE
                    # self._fastcopy(self.datadir+'punch500c.del','/dev/shm/FAL/ll.dat_{0}'.format(self.ID))
                    self.linelist = '/dev/shm/FAL/{0}/ll.dat'.format(self.ID)
                elif rlinedict['ll'] == 'master':
                    pass
                else:
                    self.linelist = '/dev/shm/FAL/{0}/ll.dat'.format(self.ID)

        # MOD ATM ->
        try:
            os.remove('/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
        except OSError:
            pass
        if self.arct_bool:
            # print('Using Arcturus model Atmosphere')
            # BOB's ARCTURUS MODEL ATM
            # self._fastcopy(self.datadir+'Arcturus_OLDpars.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
            # self._fastcopy(self.datadir+'Arcturus_NEWpars.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
            # self._fastcopy(self.datadir+'Arcturus_NEWpars_V2.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
            # self._fastcopy(self.datadir+'Arcturus_NEWpars_V3.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
            self._fastcopy(self.datadir+'Arcturus_NEWpars_V4.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
        else:
            # BOB's SOL MODEL ATM
            self._fastcopy(self.datadir+'modcaspf.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))

        # SET ALL MICROTURB VALUES TO 1.5 km/s
        # self._fastcopy(self.datadir+'modcaspf_vt15.dat','/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
        # ORIGINAL MODEL ATM
        # self._fastcopy(self.datadir+'ksolmod_syn.dat','/dev/shm/FAL/mod.dat_{0}'.format(self.ID))

        self.atmomod = '/dev/shm/FAL/{0}/mod.dat'.format(self.ID)
        
        # MOLECULES
        try:
            os.remove('/dev/shm/FAL/{0}/molecules.dat'.format(self.ID))
        except OSError:
            pass
        # BOB's MOLECULES FILE
        self._fastcopy(self.datadir+'molecules.dat','/dev/shm/FAL/{0}/molecules.dat'.format(self.ID))
        # ORIGINAL MOLECULES FILE
        # self._fastcopy(self.datadir+'molecules.dat_OLD','/dev/shm/FAL/molecules.dat_{0}'.format(self.ID))
        self.molecules = '/dev/shm/FAL/{0}/molecules.dat'.format(self.ID)
        
        # CONTINUA
        try:
            os.remove('/dev/shm/FAL/{0}/continua.dat'.format(self.ID))
        except OSError:
            pass
        # BOB's CONTINUA FILE
        self._fastcopy(self.datadir+'continua.dat','/dev/shm/FAL/{0}/continua.dat'.format(self.ID))
        # ORIGINAL CONTINUA FILE
        # self._fastcopy(self.datadir+'continua.dat_OLD','/dev/shm/FAL/continua.dat_{0}'.format(self.ID))
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
        
        # Define variables
        if self.arct_bool:
            # print('Using Arcturus model paramters')

            # self.synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
            #     "0.          0   {LINOUT:3.0f}  {TOL:7.5f}     {PRED}    00\n"
            #     "AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF  IFPREDNREAD")
            self.synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
                "0.          0   {LINOUT:3.0f}{TOL:7.5f}     {PRED}    00\n"
                "AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD")


            self.rotatevar = ("    1\n{VROT:4.3f}")

            self.macpar = "MACRO     {MACVEL:3.1f}       KM                  COMMENT FIELD"

            self.gausspar_opt = ("GAUSSIAN  130000.   RESOLUTION\n"
                             "1234567890123456789012345678901234567890")

            self.gausspar_hband = ("GAUSSIAN  130000.   RESOLUTION\n"
                             "1234567890123456789012345678901234567890")

        else:

            self.synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
                "0.          0   {LINOUT:3.0f}{TOL:7.5f}     {PRED}    72\n"
                "AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD\n"
                "  1           -0.000\n"
                "  2           -0.000\n"
                "  3           -0.000\n"
                "  4           -0.000\n"
                "  5           -0.000\n"
                "  6           -0.000\n"
                "  7           -0.000\n"
                "  8           -0.000\n"
                "  9           -0.000\n"
                " 10           -0.000\n"
                " 11           -0.000\n"
                " 12           -0.000\n"
                " 13           -0.000\n"
                " 14           -0.000\n"
                " 15           -0.000\n"
                " 16           -0.000\n"
                " 17           -0.000\n"
                " 18           -0.000\n"
                " 19           -0.000\n"
                " 20           -0.000\n"
                " 21           -0.000\n"
                " 22           -0.000\n"
                " 23           -0.000\n"
                " 24           -0.000\n"
                " 25           -0.000\n"
                " 26           -0.000\n"
                " 27           -0.000\n"
                " 28           -0.000\n"
                " 29           -0.000\n"
                " 30           -0.000\n"
                " 31           -0.000\n"
                " 32           -0.000\n"
                " 33           -0.000\n"
                " 34           -0.000\n"
                " 35           -0.000\n"
                " 36           -0.000\n"
                " 37           -0.000\n"
                " 38           -0.000\n"
                " 39           -0.000\n"
                " 40           -0.000\n"
                " 41           -0.000\n"
                " 42           -0.000\n"
                " 43           -0.000\n"
                " 44           -0.000\n"
                " 45           -0.000\n"
                " 46           -0.000\n"
                " 47           -0.000\n"
                " 48           -0.000\n"
                " 49           -0.005\n"
                " 50           -0.036\n"
                " 51           -0.173\n"
                " 52           -0.355\n"
                " 53           -0.456\n"
                " 54           -0.562\n"
                " 55           -0.679\n"
                " 56           -0.804\n"
                " 57           -0.930\n"
                " 58           -1.070\n"
                " 59           -1.152\n"
                " 60           -1.228\n"
                " 61           -1.274\n"
                " 62           -1.316\n"
                " 63           -1.356\n"
                " 64           -1.398\n"
                " 65           -1.451\n"
                " 66           -1.485\n"
                " 67           -1.519\n"
                " 68           -1.545\n"
                " 69           -1.579\n"
                " 70           -1.615\n"
                " 71           -1.662\n"
                " 72           -1.662")

            self.rotatevar = ("    1\n{VROT:4.3f}")

            self.macpar = "MACRO     {MACVEL:3.1f}       KM                  COMMENT FIELD"

            self.sincpar_opt = ("SINX/X    .047628472CM-1      .047628472   COMMENT FIELD\n"
                            "1234567890123456789012345678901234567890")

            self.gausspar_opt = ("GAUSSIAN  .071442708CM-1      .071442708   COMMENT FIELD\n"
                             "1234567890123456789012345678901234567890")

            self.sincpar_hband = ("SINX/X    .008140581CM-1      .008140581   COMMENT FIELD\n"
                            "1234567896012345678960123456789601234567890")

            self.gausspar_hband = ("GAUSSIAN  .012210871CM-1      .012210871   COMMENT FIELD\n"
                             "123456789012345678960123456789601234567890")

    def __call__(self,indict):
        """
        Function to allow a call to the class to run everything in order
        """
        # start a timer
        start_time = time.time()
        print('START CODE')
        self.synbeg(indict,verbose=False)
        synbeg_t = time.time()-start_time
        print(synbeg_t,synbeg_t)
        
        if self.rline == True:
            rlinedict = {}
            rlinedict['atoms'] = True
            rlinedict['moles'] = True
            rlinedict['TiO'] = True
            rlinedict['H2O'] = True
            rlinedict['predict'] = True
            rtype = 'rline'
        else:
            rlinedict = {}
            rtype = 'punch500'

        self.readlines(rtype=rtype,rlinedict=rlinedict,verbose=False,reuse_ll=True)
        readlines_t = time.time()-start_time
        print(readlines_t,readlines_t-synbeg_t)

        self.xnfpelsyn(verbose=False)
        xnfpelsyn_t = time.time()-start_time
        print(xnfpelsyn_t,xnfpelsyn_t-readlines_t)

        self.syn(verbose=False)
        synthe_t = time.time()-start_time
        print(synthe_t,synthe_t-xnfpelsyn_t)

        self.spectrv(verbose=False)
        spectrv_t = time.time()-start_time
        print(spectrv_t,spectrv_t-synthe_t)

        self.rotate(indict,verbose=False)
        rotate_t = time.time()-start_time
        print(rotate_t,rotate_t-spectrv_t)        

        """
        self.cleandir()
        cleandir_t = time.time()-start_time
        print(cleandir_t,cleandir_t-rotate_t)
        """

        self.broaden('ROT1',indict,broadtype='MAC')
        brdmac_t = time.time()-start_time
        # print(brdmac_t,brdmac_t-cleandir_t)

        self.broaden('ROT1_mac',indict,broadtype='INSTRUMENT')
        brdins_t = time.time()-start_time
        print(brdins_t,brdins_t-brdmac_t)

        print('END CODE')       

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

    def synbeg(self,indict,verbose=None):
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
        # read in information from input dictionary
        WSTART = indict['WSTART']
        WEND = indict['WEND']
        RESOL = indict['RESOL']
        PRED = indict['PRED']
        if 'LINOUT' in indict.keys():
            LINOUT = indict['LINOUT']
        else:
            LINOUT = 10
        TOL = 1e-3
        # TOL = 0.0
        
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
        synbegvar_i = self.synbegvar.format(WSTART=WSTART,WEND=WEND,RESOL=RESOL,PRED=PRED,TOL=TOL,LINOUT=LINOUT)
        self.synbegout = self._callpro("synbeg",synbegvar_i,verbose=verbose)
        if verbose:
            print("Finished SynBeg")

        # make copy of files in INT directory so that everything can be restarted fresh
        if not os.path.exists('/dev/shm/FAL/{0}/INT'.format(self.ID)):
            os.makedirs('/dev/shm/FAL/{0}/INT'.format(self.ID))
        fortlist = glob.glob('fort.*')
        for ft in fortlist:
            self._fastcopy('/dev/shm/FAL/{0}/{1}'.format(self.ID,ft),
                           '/dev/shm/FAL/{0}/INT/{1}'.format(self.ID,ft))

        return (self.synbegout,self.ID)

    def readlines(self,rtype=None,rlinedict={},verbose=None,reuse_ll=True):
        """
        function that determines if you should use a punch500 or 
        read in all the individual lines
        
        rtype = None defaults to using the punch500_new.del file
              = 'punch500' -> punch500_new.del
              = 'rline' -> read in all of the lines
        """

        if (rtype == None or rtype == "punch500"):
            self.rgfalldel(verbose=verbose,reuse_ll=reuse_ll)
        elif (rtype == 'master'):
            self.rmaster(rlinedict['MASTERLL'],verbose=verbose)
        elif (rtype == 'rline'):
            self.rlinefunc(rlinedict,verbose=verbose)
        else:
            raise IOError('Did not understand read line type rtype={0}'.format(rtype))

    def rmaster(self,MASTERLL=None,MASTERMOLLL=None,verbose=None):
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
        # link master line list file
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)

        # Master lines file is gigantic, so don't copy into memory just sym link it
        if MASTERLL == None:
            MASTERLL = '/work/02349/cargilpa/FAL/MASTERLL/HBAND/CargileLL_1400_1900.bin'

        os.symlink(MASTERLL,'fort.11')

        # run read master list program
        self.rmasterout = self._callpro("rpunchbin",verbose=verbose)

        # link master line list file
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)

        if MASTERMOLLL == None:
            MASTERMOLLL = '/work/02349/cargilpa/FAL/MASTERLL/HBAND/H2OLL_1400_1900.bin'

        # H2O+TiO lines file is gigantic, so don't copy into memory just sym link it
        os.symlink(MASTERMOLLL,'fort.11')

        # run read master list program
        self.rmasterout = self._callpro("rpunchbin",verbose=verbose)


        # # check to make sure all input/output files are right
        # filesdict = {}
        # filesdict['infiles'] = {'fort.11':'bin'}
        # filesdict['writeto'] = ({'fort.12':'bin','fort.14':'bin','fort.19':'bin',
        #                           'fort.20':'bin','fort.93':'bin'})
        # filesdict['newfiles'] = {}
        # try:
        #     assert self._fileprep(filesdict) == None
        #     pass
        # except AssertionError:
        #     raise IOError("Something wrong with Input/Output files")


        # move fort.14 to fort.99 and clean up
        self._mvsym('fort.14','fort.99')
        
        return (self.rmasterout,self.ID)


    def rgfalldel(self,verbose=None,reuse_ll=False):
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
            self._makesym(self.linelist,'fort.11')
        
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
                # molefiles = ([
                #         'alopatrascu.asc', # AlO
                #         'nah.dat', # NaH
                #         'voax.asc','vobx.asc','vocx.asc', # VO
                #         'fehfx.dat', # FeH
                #         'h2bx.dat','h2cx.dat','h2xx.dat', # H2
                #         'hdxx.dat', #HD
                #         'mgh.dat', #MgH
                #         'mgodaily.dat', #MgO
                #         'nhax.dat','nhca.dat', # NH
                #         'chjorg.dat', #CH
                #         'cnaxbrooke.dat','cnbxbrooke.dat','cnxx12brooke.dat', # CN
                #         'c2ax.dat','c2ba.dat','c2da.dat','c2ea.dat', #C2
                #         'coax.dat','coxx.dat', # CO
                #         'ohax.dat','ohxx.dat', # OH
                #         'sihax.dat', # SiH
                #         'sioax.dat','sioex.dat','sioxx.dat', #SiO
                #         'crhax.dat', # CrH
                #         'cah.dat', # CaH
                #     ])
                # read molecular files
                self.rmolecascout = {}
                [self.readmol(mf,verbose=verbose) for mf in self.molefiles]
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
        # write TiO line list files into fort.11 and fort.48
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        if self.savetomem:
            self._makesym('/dev/shm/FAL/{0}/schwenke.bin'.format(self.ID),'fort.11')
        else:
            os.symlink(self.moldatadir+'schwenke.bin','fort.11')
            os.symlink(self.moldatadir+'schwenke.bin','/dev/shm/FAL/{0}/fort.11'.format(self.ID))        
        if os.path.isfile("fort.48"):
            self._rmsym('fort.48',verbose=verbose)
        if self.savetomem:
            self._makesym('/dev/shm/FAL/{0}/eschwenke.bin'.format(self.ID),'fort.48')
        else:
            os.symlink(self.moldatadir+'eschwenke.bin','fort.48')
            os.symlink(self.moldatadir+'eschwenke.bin','/dev/shm/FAL/{0}/fort.48'.format(self.ID))        

        # print("Running RSchwenk on TiO")
        self.rmoleout_tio = self._callpro("rschwenk",verbose=verbose)
        # print("Finished RSchwenk")

        if os.path.isfile("fort.48"):
            os.unlink('fort.48')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.48'.format(self.ID)):
            os.remove('/dev/shm/FAL/{0}/fort.48'.format(self.ID))
    
    def readmol_H2O(self,verbose=None):
        """
        helper function to do molecule H2O file read in
        """
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        if self.savetomem:
            self._makesym('/dev/shm/FAL/{0}/h2oslowfix.bin'.format(self.ID),'fort.11')
        else:
            os.symlink(self.moldatadir+'h2oslowfix.bin','fort.11')
            os.symlink(self.moldatadir+'h2oslowfix.bin','/dev/shm/FAL/{0}/fort.11'.format(self.ID))        
        if os.path.isfile("fort.48"):
            self._rmsym('fort.48',verbose=verbose)
        if self.savetomem:
            self._makesym('/dev/shm/FAL/{0}/eh2opartridge.bin'.format(self.ID),'fort.48')
        else:
            os.symlink(self.moldatadir+'eh2opartridge.bin','fort.48')
            os.symlink(self.moldatadir+'eh2opartridge.bin','/dev/shm/FAL/{0}/fort.48'.format(self.ID))        
        # print("Running RH2OFast on H2O")
        self.rmoleout_h2o = self._callpro("rh2oslow",verbose=verbose)
        # print("Finished RH2OFast")

        if os.path.isfile("fort.48"):
            os.unlink('fort.48')
        if os.path.isfile('/dev/shm/FAL/{0}/fort.48'.format(self.ID)):
            os.remove('/dev/shm/FAL/{0}/fort.48'.format(self.ID))

    def readmol(self,molfile,verbose=None):
        """
        helper function to do molecule file read in
        """
        # write molfile line list into fort.11
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        if self.savetomem:
            self._makesym('/dev/shm/FAL/{0}/{1}'.format(self.ID,molfile),'fort.11')
        else:
            os.symlink(self.moldatadir+'{0}'.format(molfile),'fort.11')
            os.symlink(self.moldatadir+'{0}'.format(molfile),'/dev/shm/FAL/{0}/fort.11'.format(self.ID))        
        
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
        # os.symlink('/work/02349/cargilpa/FAL/DATA/gfpred29dec2014.bin','fort.11')
        os.symlink('/work/02349/cargilpa/FAL/DATA/gfpred27sep15.bin','fort.11')
        # print("Running RPredict")
        self.rmoleout = self._callpro("rpredict",verbose=verbose)
        # print("Finished RPredict")

    def ratomic(self,verbose=None):
        """
        helper function to read in atomic lines
        """
        # write atomic line list into fort.11
        if os.path.isfile("fort.11"):
            self._rmsym('fort.11',verbose=verbose)
        if self.savetomem:
            self._makesym(self.atomiclinelist,'fort.11')
        else:            
            os.symlink('/work/02349/cargilpa/FAL/DATA/gfall25jul15.dat','fort.11')
            os.symlink('/work/02349/cargilpa/FAL/DATA/gfall25jul15.dat','/dev/shm/FAL/{0}/fort.11'.format(self.ID))

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

    def rotate(self,indict,verbose=None):
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
        rotatestr = self.rotatevar.format(VROT=indict['VROT'])
        self.rotateout = self._callpro("rotate",rotatestr,verbose=verbose)
        # print("Finished Rotate")

        return (self.rotateout,self.ID)
    
    def broaden(self,inspec,indict,broadtype=None,WLreg=None,write=False,verbose=None):
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
            self.broadout = self._callpro(
                "broadenx",inputstr=self.macpar.format(MACVEL=indict['MACVEL']),
                verbose=verbose)

        elif broadtype=="INSTRUMENT":
            tag = 'inst'
            # determine if arcturus or not
            if self.arct_bool:
                # determine which instrument settings should be used
                if WLreg=='OPT':
                    gausspar = self.gausspar_opt
                elif WLreg=='HBAND':
                    gausspar = self.gausspar_hband
                else:
                    gausspar = self.gausspar_opt

                self.broadout2 = self._callpro("broadenx",inputstr=gausspar,verbose=verbose)
            else:
                # determine which instrument settings should be used
                if WLreg=='OPT':
                    sincpar = self.sincpar_opt
                    gausspar = self.gausspar_opt
                elif WLreg=='HBAND':
                    sincpar = self.sincpar_hband
                    gausspar = self.gausspar_hband
                else:
                    sincpar = self.sincpar_opt
                    gausspar = self.gausspar_opt

                # first run sinc
                self.broadout1 = self._callpro("broadenx",inputstr=sincpar,verbose=verbose)
                # mv fort.22 into fort.21 as into into gauss
                os.unlink('fort.21')
                os.remove('/dev/shm/FAL/{0}/fort.21'.format(self.ID))
                self._mvsym('fort.22','fort.21')
                ff = open('/dev/shm/FAL/{0}/fort.22'.format(self.ID),'wb')
                f = np.fromfile(ff,count=0)
                ff.close()
                os.symlink('/dev/shm/FAL/{0}/fort.22'.format(self.ID),'fort.22')
                self.broadout2 = self._callpro("broadenx",inputstr=gausspar,verbose=verbose)
        else:
            raise ValueError("Did not understand broadening type (MAC/INSTRUMENT)")
        
        # move final fort.22 into an output file with the right tag
        if os.path.exists(inspec+'_{0}'.format(tag)):
            os.remove(inspec+'_{0}'.format(tag))
        self._mvsym('fort.22',inspec+'_{0}'.format(tag))

    def transmit(self):
        pass

    def writespec(self,inspec,local=False,verbose=None):
        """
        Function to write out binary spectrum to ascii
        """
        
        # mv input file into fort.21
        try:
            os.remove('/dev/shm/FAL/fort.1_'+self.ID)
        except OSError:
            pass
        try:
            os.remove('fort.1_'+self.ID)
        except OSError:
            pass

        try:
            os.remove('/dev/shm/FAL/fort.2_'+self.ID)
        except OSError:
            pass
        try:
            os.remove('fort.2_'+self.ID)
        except OSError:
            pass

        try:
            os.remove('/dev/shm/FAL/fort.3_'+self.ID)
        except OSError:
            pass
        try:
            os.remove('fort.3_'+self.ID)
        except OSError:
            pass

        try:
            os.remove('/dev/shm/FAL/fort.4_'+self.ID)
        except OSError:
            pass
        try:
            os.remove('fort.4_'+self.ID)
        except OSError:
            pass

        self._makesym('/dev/shm/FAL/'+inspec+'_'+self.ID,'fort.1_'+self.ID)

        # check to make sure all input/output files are right
        filesdict = {}
        filesdict['infiles'] = {'fort.1':'bin'}
        filesdict['writeto'] = {}
        filesdict['newfiles'] = {'fort.2':'ascii','fort.3':'ascii','fort.4':'ascii'}

        try:
            assert self._fileprep(filesdict) == None
            pass
        except AssertionError:
            raise IOError("Something wrong with Input/Output files")

        print('Writing Spectrum {0} ...'.format(inspec+'_'+self.ID))
        self.writespecout = self._callpro('syntoascanga',verbose=verbose)
        self._mvsym('fort.2_'+self.ID,'spec_'+inspec+'_'+self.ID+'.asc') 
        self._mvsym('fort.3_'+self.ID,'line_'+inspec+'_'+self.ID+'.asc')
        self._mvsym('fort.4_'+self.ID,'head_'+inspec+'_'+self.ID+'.asc')

        if local:
            subprocess.check_call(
                'mv /dev/shm/FAL/spec_'+inspec+'_'+self.ID+'.asc ./spec_'+inspec+'_'+self.ID+'.asc',shell=True)
            subprocess.check_call(
                'mv /dev/shm/FAL/head_'+inspec+'_'+self.ID+'.asc ./head_'+inspec+'_'+self.ID+'.asc',shell=True)
            subprocess.check_call(
                'mv /dev/shm/FAL/line_'+inspec+'_'+self.ID+'.asc ./line_'+inspec+'_'+self.ID+'.asc',shell=True)
        
    def at12tosyn(self):
        pass

    def reset(self):
        """
        RESET BUTTON
        """
        # remove all all old files except fort.11
        fortlist = glob.glob('fort.*') + glob.glob('ROT*')
        if 'fort.11' in fortlist:
            fortlist.remove('fort.11')

        for ft in fortlist:
            self._rmsym(ft)

        ## remove any ascii output files
        # outfiles = ['headinfo.dat','lineinfo.dat','specfile.dat']
        # for ofile in outfiles:
        #   os.remove(ofile)

        # copy files for INT back into memory and link them
        fortlist = ['fort.10','fort.12','fort.14','fort.19','fort.20','fort.93']
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
