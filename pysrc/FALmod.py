from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from astropy.table import Table
import FALsynthepkg
import FALGlue
import shutil, os, subprocess
import numpy as np
import sys,glob,time,shutil

class FALmod(object):
    """
    Class to take delta line parameters and return a synthesized
    spectrum from SYNTHE to use in python
    """    
    def __init__(self,ID=None,starpars='Sun',waverange=None):

        # set ID
        if ID != None:
            self.IDraw = ID
            self.ID = str(ID).rjust(8,str(0))
        else:
            import uuid
            self.ID = uuid.uuid4().hex[:8]
            self.IDraw = None

        # set starpars if starpars='Sun' or starpars='Arcturus'
        if starpars == 'Sun':
            self.starpars = {}
            self.starpars['VROT'] = -2.02
            self.starpars['MACVEL'] = 1.5
            self.starpars['OBJECT'] = 'Sun'
        elif starpars == 'Arcturus':
            self.starpars = {}
            self.starpars['VROT'] = 1.0
            self.starpars['MACVEL'] = 3.0
            self.starpars['OBJECT'] = 'Arcturus'
        else:
            self.starpars = starpars

        if waverange == None:
            self.starpars['WSTART'] = 600.0
            self.starpars['WEND'] = 601.0
        else:
            self.starpars['WSTART'] = waverange[0]
            self.starpars['WEND'] = waverange[1]

        # create the glue  
        self.glue = FALGlue.glue()

        # initialize SYNTHE package
        self.SYNTHE = FALsynthepkg.synthe(ID=ID,verbose=verbose,clobber=True,starpars=self.starpars)

        # make dir in memory if not already there
        if not os.path.exists('/dev/shm/FAL/{0}'.format(self.ID)):
            os.makedirs('/dev/shm/FAL/{0}'.format(self.ID))

        # make local dir to work in
        if not os.path.exists('{0}'.format(self.ID)):
            os.makedirs('{0}'.format(self.ID))

        # define parent directory
        self.parentdir = os.getcwd()

    def runsynthe(self,**kwargs):
        """
        Call SYNTHE
        """
        # sort out keywords
        # timing information printed out
        self.timeit = kwargs.get('timeit',False)
        # type of line list reader
        linelist = kwargs.get('linelist','master')
        # set any offsets in line list
        parr = kwargs.get('parr',None)
        # copy synbeg files into INT/ so that they can be reused
        archive = kwargs.get('archive',False)
        # multiply by transmission spec given as input as transspec
        transspec = kwargs.get('transspec',False)
        # write spectrum, linelist, and header out as ascii
        writespec = kwargs.get('writespec',False)
        # do work verbosely
        verbose = kwargs.get('verbose',False)


        # change into working directory
        os.chdir('{0}'.format(self.ID))

        # start if timer if needed
        if self.timeit:
            self.starttime = time.time()
            self.lasttime = self.starttime

        # call synbeg if readline != readold
        if (verbose == True or verbose == 'synbeg'):
            verbose_i = True
        else:
            verbose_i = False

        if linelist == 'readlast':
            if not os.path.exists('/dev/shm/FAL/{0}/INT'.format(self.ID)):
                raise IOError("Pro: {1} --> WARNING: COULD NOT FIND INITIAL FILES!!!! {0:7.5f} s".format(time.time()-self.starttime,self.IDraw))
            else:
                self.SYNTHE.reset()
        else:
            self.SYNTHE.synbeg(self.starpars,verbose=verbose_i)

        # read in line list
        if (verbose == True or verbose == 'readlines'):
            verbose_i = True
        else:
            verbose_i = False
        ll = self._readline(linelist,verbose_i)

        # adjust any delta values if parr has been provided
        if (verbose == True or verbose == 'adjustlines'):
            verbose_i = True
        else:
            verbose_i = False
        if parr != None:
            self._adjustpar(parr,verbose_i)

        # do synthesis calc
        if (verbose == True) or ('synthesis' in verbose):
            verbose_i = verbose
        else:
            verbose_i = False
        self._synthesis(verbose_i)

        # do broadening
        if (verbose == True) or ('broaden' in verbose):
            verbose_i = verbose
        else:
            verbose_i = False
        binspecname = self._broaden(verbose_i)

        # do specout to get final spectrum
        if (verbose == True or verbose == 'specout'):
            verbose_i = True
        else:
            verbose_i = False
        outspec,newll = self._specout(binspecname,verbose_i)

        if archive:
            # write newll into INT/ for archiving purposes
            if os.path.exists('fort.11'):
                os.unlink('fort.11')
            if os.path.exists('/dev/shm/FAL/{0}/fort.11'.format(self.ID)):
                os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))
            newll_st = self.glue.con_nptolp(newll)
            self.glue.writelp(newll_st,'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
            os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')
            self.SYNTHE.archive()

        if transspec != False:
            # do multiply transmission spectrum
            outspec = self._multiplytrans(outspec,transspec)

        if writespec:
            # write out spec, line list, and header to ascii
            self._write(binspecname)

        # cd back into parent directory
        os.chdir(self.parentdir)

        if timeit:
            print("Pro: {1} --> FINISHED FALmod -- Full time: {0:7.5f} s".format(time.time()-self.starttime,self.IDraw))

        return (outspec,newll)

    def __del__(self):
        os.chdir(self.parentdir)

    def _readline(self,linelist,verbose_i=False):
        # read the line list appropriate for the user call

        if linelist == 'readall':
            # read all individual line lists
            rlinedict = {"atoms":True,"moles":True,"H2O":True,"TiO":True} # atoms, molecules + H2O & TiO
            self.SYNTHE.readlines(rtype='readall',rlinedict=rlinedict,verbose=verbose_i)

            # # remove the various line files in memory to save resources
            # datfiles = glob.glob('/dev/shm/FAL/{0}/*.dat'.format(self.ID))
            # binfiles = glob.glob('/dev/shm/FAL/{0}/*.bin'.format(self.ID))
            # savefiles = ['he1tables.dat','continua.dat','molecules.dat','mod.dat']
            # addfiles  = ['voax.asc','vobx.asc','vocx.asc']
            # [datfiles.remove('/dev/shm/FAL/{0}/{1}'.format(self.ID,sf)) for sf in savefiles]
            # [datfiles.append('/dev/shm/FAL/{0}/{1}'.format(self.ID,af)) for af in addfiles]
            # [os.remove(dfil) for dfil in datafiles if os.path.isfile(dfil)]
            # [os.remove(bfil) for bfil in binfiles if os.path.isfile(bfil)]

            if self.timeit:
                print("Pro: {1} --> Read in all line lists -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()

        elif linelist == 'readmaster':
            # read the masterline lists (cargile or kurucz plus H2O+TiO)
            self.SYNTHE.readlines(rtype='readmaster',verbose=verbose_i)
            if self.timeit:
                print("Pro: {1} --> Read in master line list -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()

        elif linelist == 'readlast':
            # read the previously created line list in directory
            self.SYNTHE.readlines(rtype='readlast',verbose=verbose_i)
            if self.timeit:
                print("Pro: {1} --> Read in lines -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()

        else:
            # linelist equal to a path to user defined line list 
            self.SYNTHE.readlines(rtype=linelist,verbose=verbose_i)
            if self.timeit:
                print("Pro: {1} --> Read in user defined line list {2} -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw,linelist))
                self.lasttime = time.time()


    def _adjustpar(self,parr_i,ll=None,verbose=False):
        # check to see if parr is one table of parameters (for the whole line list)
        # or if it is a dictionary of parameters and line indices
        if isinstance(parr,dict):
            parr = parr_i['parr']
            lineind = parr_i['lineind']
        else:
            parr = parr_i

        # check to see if line list has been defined in self
        if ll == None:
            ll_i = self.glue.readlp_raw('/dev/shm/FAL/{0}/fort.11'.format(self.ID))
            ll = self.glue.con_lptonp_raw(ll_i)

        # copy orginal linelist into a working list
        self.ll = ll.copy()

        # delete old fort.11 file
        os.unlink('fort.11')
        os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))

        # define free parameters and stick them into line list
        # FREE POSIBLE PARS: DWL, DGFLOG, and one or all of DGAMMAR, DGAMMAS, DGAMMAW

        # the case where all lines are free (i.e., no lineind index array)
        if lineind==None:
            del self.ll['DWL','DGFLOG','DGAMMAR','DGAMMAS','DGAMMAW']
            self.ll.add_columns(p.columns.values())
        else:
            # else just replace the parameters for the lines that are indictated in lineind
            for pind,llind in enumerate(lineind):
                self.ll['DWL'][llind] = p['DWL'][pind]
                self.ll['DGFLOG'][llind] = p['DGFLOG'][pind]
                self.ll['DGAMMAR'][llind] = p['DGAMMAR'][pind]
                self.ll['DGAMMAS'][llind] = p['DGAMMAS'][pind]
                self.ll['DGAMMAW'][llind] = p['DGAMMAW'][pind]

        # convert the table into the correct string format
        lpfmttab = self.glue.con_nptolp(self.ll)

        # write out the ascii line list
        self.glue.writelp(lpfmttab,'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
        os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')

    def _synthesis(self,verbose=False):
        # now run the synthe steps

        # -- run xnfpelsyn --
        if (verbose == True or verbose == 'synthesis:xnfpelsyn'):
            verbose_i = True
        else:
            verbose_i = False
        self.SYNTHE.xnfpelsyn(verbose=verbose_i)
        if self.timeit:
            print("Pro: {1} --> XNFPELSYN -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

        # -- do synthe --
        if (verbose == True or verbose == 'synthesis:synthe'):
            verbose_i = True
        else:
            verbose_i = False
        self.SYNTHE.syn(verbose=verbose_i,speed=self.speed)
        if self.timeit:
            print("Pro: {1} --> SYNTHE -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

        # -- do spectrv --
        if (verbose == True or verbose == 'synthesis:spectrv'):
            verbose_i = True
        else:
            verbose_i = False
        self.SYNTHE.spectrv(verbose=verbose_i,tau=self.tau_i)
        if self.timeit:
            print("Pro: {1} --> SPECTRV -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

    def _broaden(self,verbose=False):
        # -- do rotate --
        if (verbose == True or verbose == 'broaden:rotate'):
            verbose_i = True
        else:
            verbose_i = False
        self.SYNTHE.rotate(self.starpars['VROT'],verbose=verbose_i)
        if self.timeit:
            print("Pro: {1} --> ROTATE -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

        # -- check if the user wants broadening --
        if self.starpars['BROADEN'] == -1:
            # no macrovel applied, just copy rotated spectrum to output
            self.SYNTHE._makesym('/dev/shm/FAL/{0}/{1}'.format(self.ID,'ROT1'),'ROT1_mac_inst')
            print("Pro: {1} --> No Broadening Applied -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()
    
        else:
            # -- do broaden for macroturblence --
            if (verbose == True or verbose == 'broaden:broad_mac'):
                verbose_i = True
            else:
                verbose_i = False
            self.SYNTHE.broaden('ROT1',self.starpars['MACVEL'],broadtype='MAC',verbose=verbose_i)
            if self.timeit:
                print("Pro: {1} --> BROADEN MACRO -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()

            # -- do broaden for instrumental --
            if (verbose == True or verbose == 'broaden:broad_inst'):
                verbose_i = True
            else:
                verbose_i = False
            self.SYNTHE.broaden('ROT1_mac',broadtype='INSTRUMENT',verbose=verbose_i)
            if self.timeit:
                print("Pro: {1} --> BROADEN INSTR -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()

    def _specout(self,infile):
        # read in binary output spectrum

        # pull run info from RUNINFO.dat file
        runinfo = {}        
        with open('RUNINFO.dat','r') as runinfoff:
            runinfoarr = runinfoff.readlines()
            runinfo_i = runinfoarr[1].split()
            runinfo['LENGTH'] = int(runinfo_i[0])
            runinfo['RATIO'] = float(runinfo_i[1])
            runinfo['WBEGIN'] = float(runinfo_i[2])
            runinfo['DWLBEG'] = float(runinfo_i[3])
            runinfo['WLLAST'] = float(runinfo_i[4])
            runinfo['DWLLAST'] = float(runinfo_i[5])

        # pull line info from LINEINFO.dat file
        lineinfo = {}
        with open('LINEINFO.dat','r') as lineinfoff:
            lineinfoarr = lineinfoff.readlines()
            lineinfo_i = lineinfoarr[1]
            lineinfo['NLINES'] = int(lineinfo_i)

        outspec,newll = self.glue.readspecbin(infile,NWL=runinfo['LENGTH'],NLINES=lineinfo['NLINES'])
        # turn these dictionaries into Table for ease of use
        outspec = Table(outspec)
        newll = Table(newll)

        if self.timeit:
            print("Pro: {1} --> Read in binary spectrum -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()
        return (outspec,newll)


    def _multiplytrans(self,inspec,trans):
        # multiply inspec you input transmissino spectrum, should be on same wavelength scale (may change this to be more general)
        outspec = np.multiply(inspec['QMU1'],trans['FLUX'])
        if self.timeit:
            print("Pro: {1} --> Multiply Transmission Spectrum -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()
        return outspec

    def _writeout(self,inbinary):
        # write out ascii version of spectrum, line list, and header file

        self.SYNTHE.writespec(inbinary)
        if self.timeit:
            print("Pro: {1} --> Write ASCII files -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()



# class FALmod(object):
#     """
#     Class to take delta line parameters and return a synthesized
#     spectrum from SYNTHE to use in python
#     """    
#     def __init__(self,startll=None,ID=None,indict=None,rlinestop=False,injectlines=None,arct_bool=False,linelistcheck=False):

#         if ID != None:
#             self.IDraw = ID
#             self.ID = str(ID).rjust(8,str(0))
#         else:
#             import uuid
#             self.ID = uuid.uuid4().hex[:8]
#             self.IDraw = None

#         # set arcturus boolean
#         self.arct_bool = arct_bool

#         # create the glue        
#         self.glue = FALGlue.glue()

#         # make dir in memory if not already there
#         if not os.path.exists('/dev/shm/FAL/{0}'.format(self.ID)):
#             os.makedirs('/dev/shm/FAL/{0}'.format(self.ID))

#         # make local dir to work in
#         if not os.path.exists('{0}'.format(self.ID)):
#             os.makedirs('{0}'.format(self.ID))

#         # define parent directory
#         self.parentdir = os.getcwd()

#         # change into working directory
#         os.chdir('{0}'.format(self.ID))

#         # set up boolean flags for SYNTHE
#         if (startll == None) or (startll == 'rline'):
#             self.LLbool = False
#             self.rlinebool = True
#             # self.rlinedict = {"atoms":True} # just atoms
#             # self.rlinedict = {"atoms":True,"moles":True} # atoms and molecules
#             # self.rlinedict = {"atoms":True,"moles":True,"TiO":True} # no water
#             self.rlinedict = {"atoms":True,"moles":True,"H2O":True,"TiO":True} # atoms, molecules + H2O & TiO
#             # self.rlinedict = {"atoms":True,"H2O":True} # atoms + H2O

#             if ("PRED" in indict.keys()):
#                 if indict["PRED"] == True:
#                     self.rlinedict['predict'] = True
#             if (injectlines != None):
#                 self.rlinedict['injectlines'] = injectlines

#         elif (startll == 'punch500'):
#             self.LLbool = False
#             self.rlinebool = False
#             self.rlinedict = {'ll':'punch500'}
#         elif (startll == 'master'):
#             self.LLbool = False
#             self.rlinebool = False
#             self.rlinedict = {'ll':'master'}
#         else:
#             self.LLbool = True
#             self.rlinebool = False
#             self.rlinedict = {}

#         # read in a file that indicates which of the gamma parameters is dominate
#         # default (for lines that are not in self.ll) = DGAMMAW
#         gammafile = Table.read('/work/02349/cargilpa/FAL/DATA/MASTERGAMMASWITCH_V3.fits',format='fits')

#         # initialize SYNTHE Python Wrapper
#         SYNTHE = synthepkg.synthe(ID=ID,verbose=False,clobber=True,rline=self.rlinebool,LL=self.LLbool,rlinedict=self.rlinedict,arct_bool=self.arct_bool)

#         # build indict if not stated
#         if indict == None:
#             self.indict = {'WSTART':517.6,'WEND':517.8,'RESOL':1000000.0,'VROT':-2.02,'MACVEL':1.5,'PRED':0}
#         else:
#             self.indict = indict
#             SYNTHE.synbeg(self.indict,verbose=False)

#         # determine the WLreg settings
#         if (indict['WSTART'] > 472.0) & (indict['WEND'] < 1000.00):
#             self.WLreg = 'OPT'
#         elif (indict['WSTART'] > 1000.00) & (indict['WEND'] < 3000.00):
#             self.WLreg = 'HBAND'
#         else:
#             self.WLreg = 'OPT'
#             print("Pro: {0} --> WARNING: using default OPTICAL instrument settings!".format(self.IDraw))

#         # now figure out how to deal with line list and write into self
#         if (startll == None) or (startll == 'rline'):
#             print("Pro: {0} --> Running Full Set of Line List Files...".format(self.IDraw))
            
#             starttime = time.time()
#             makeverb = False
#             SYNTHE.readlines(rtype='rline',rlinedict=self.rlinedict,verbose=False,reuse_ll=False)
#             # remove the various line files in memory to save resources
#             datfiles = glob.glob('/dev/shm/FAL/{0}/*.dat'.format(self.ID))
#             datfiles.remove('/dev/shm/FAL/{0}/he1tables.dat'.format(self.ID))
#             datfiles.remove('/dev/shm/FAL/{0}/continua.dat'.format(self.ID))
#             datfiles.remove('/dev/shm/FAL/{0}/mod.dat'.format(self.ID))
#             datfiles.remove('/dev/shm/FAL/{0}/molecules.dat'.format(self.ID))
#             datfiles.append('/dev/shm/FAL/{0}/voax.asc'.format(self.ID))
#             datfiles.append('/dev/shm/FAL/{0}/vobx.asc'.format(self.ID))
#             datfiles.append('/dev/shm/FAL/{0}/vocx.asc'.format(self.ID))
            
#             for dfil in datfiles:
#                 if os.path.isfile(dfil):
#                     os.remove(dfil)

#             binfiles = glob.glob('/dev/shm/FAL/{0}/*.bin'.format(self.ID))
#             for bfil in binfiles:
#                 if os.path.isfile(bfil):
#                     os.remove(bfil)

#             SYNTHE.xnfpelsyn(verbose=False)
#             SYNTHE.syn(verbose=False,speed='slow')
#             SYNTHE.spectrv(verbose=False)
#             SYNTHE.rotate(self.indict,verbose=False)

#             if os.path.isfile("fort.1"):
#                 os.unlink('fort.1')
#             shutil.copy('/dev/shm/FAL/{0}/ROT1'.format(self.ID),'fort.1')
#             subprocess.check_call('/home1/02349/cargilpa/FAL/PYTHON/bin/syntoascanga.exe',
#                     shell=True,stdout=open(os.devnull,"w"))
#             ll_i = self.glue.readlp_raw('lineinfo.dat')
#             ll = self.glue.con_lptonp_raw(ll_i)


#             # SYNTHE.broaden('ROT1',self.indict,broadtype='MAC',verbose=False)
#             # SYNTHE.broaden('ROT1_mac',self.indict,broadtype='INSTRUMENT',WLreg=self.WLreg,verbose=False)
#             # shutil.copy('/dev/shm/FAL/{0}/ROT1_mac_inst'.format(self.ID),'fort.1')
#             # subprocess.check_call('/home1/02349/cargilpa/FAL/PYTHON/bin/syntoascanga.exe',
#             #     shell=True,stdout=open(os.devnull,"w"))
#             # ll_i = self.glue.readlp_raw('lineinfo.dat')
#             # ll = self.glue.con_lptonp_raw(ll_i)

#             # copy orginal linelist into a working list
#             ll.sort('WL')
#             self.ll = ll.copy()

#             # when using full set of line list, change ll.dat into fort.11 just for coding reasons

#             # delete old fort.11 file
#             if os.path.isfile("fort.11"):
#                 os.unlink('fort.11')
#             if os.path.isfile('/dev/shm/FAL/{0}/fort.11'.format(self.ID)):
#                 os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))

#             # convert the table into the correct string format
#             lpfmttab = self.glue.con_nptolp(self.ll)

#             # copy line list into ll.dat file
#             self.glue.writelp(lpfmttab,'/dev/shm/FAL/{0}/ll.dat'.format(self.ID))

#             # write out the ascii line list
#             self.glue.writelp('/dev/shm/FAL/{0}/ll.dat'.format(self.ID),'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
#             os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')

#             # reinitialize SYNTHE Python Wrapper as long as rlinestop isn't triggered
#             if rlinestop==False:
#                 SYNTHE.reset()
#                 self.SYNTHE = synthepkg.synthe(ID=ID,verbose=False,clobber=True,rline=False,LL=True,rlinedict={},arct_bool=self.arct_bool)
#                 self.indict['LINEOUT'] = 30
#                 self.SYNTHE.synbeg(self.indict,verbose=False)

#             print("Pro: {0} --> Full set of lines took {1}s".format(self.IDraw,time.time()-starttime))

#         elif startll == 'punch500':
#             print("Pro: {0} --> Running Default Punch500 File...".format(self.IDraw))
#             makeverb = False
#             SYNTHE.readlines(rtype='punch500',rlinedict=self.rlinedict,verbose=False,reuse_ll=False)
#             SYNTHE.xnfpelsyn(verbose=False)
#             SYNTHE.syn(verbose=False,speed='fast')
#             SYNTHE.spectrv(verbose=False)
#             SYNTHE.rotate(self.indict,verbose=False)
#             os.unlink('fort.1')
#             shutil.copy('/dev/shm/FAL/{0}/ROT1'.format(self.ID),'fort.1')
#             subprocess.check_call('/home1/02349/cargilpa/FAL/PYTHON/bin/syntoascanga.exe',
#                     shell=True,stdout=open(os.devnull,"w"))
#             ll_i = self.glue.readlp_raw('lineinfo.dat')
#             ll = self.glue.con_lptonp_raw(ll_i)

#             # copy orginal linelist into a working list
#             self.ll = ll.copy()

#             # copy line list into ll.dat file
#             shutil.copy('lineinfo.dat','/dev/shm/FAL/{0}/ll.dat'.format(self.ID))

#             # reinitialize SYNTHE Python Wrapper
#             SYNTHE.reset()
#             self.SYNTHE = synthepkg.synthe(ID=ID,verbose=False,clobber=True,rline=False,LL=True,rlinedict={},arct_bool=self.arct_bool)
#             self.SYNTHE.synbeg(self.indict,verbose=False)

#         elif startll == 'master':
#             self.rlinedict['MASTERLL'] = '/work/02349/cargilpa/FAL/MASTERLL/HBAND/CargileLL_1400_1900.bin'
#             # self.rlinedict['MASTERLL'] = '/work/02349/cargilpa/FAL/MASTERLL/HBAND/KuruczLL_1400_1900.bin'
#             # self.rlinedict['MASTERLL'] = '/work/02349/cargilpa/FAL/SIMPLEKRZ/new_kurucz/data/FullLL_1400_1900.bin'
#             self.rlinedict['MASTERMOLLL'] = '/work/02349/cargilpa/FAL/MASTERLL/HBAND/KuruczH2OLL_1400_1900.bin'
#             print("Pro: {0} --> Running Master Line List File: {1} & {2}...".format(self.IDraw,self.rlinedict['MASTERLL'],self.rlinedict['MASTERMOLLL']))
#             starttime = time.time()
#             makeverb = False
#             SYNTHE.readlines(rtype='master',rlinedict=self.rlinedict,verbose=False,reuse_ll=False)
#             SYNTHE.xnfpelsyn(verbose=False)
#             SYNTHE.syn(verbose=False,speed='slow')
#             SYNTHE.spectrv(verbose=False)
#             SYNTHE.rotate(self.indict,verbose=False)
#             os.unlink('fort.1')
#             shutil.copy('/dev/shm/FAL/{0}/ROT1'.format(self.ID),'fort.1')
#             subprocess.check_call('/home1/02349/cargilpa/FAL/PYTHON/bin/syntoascanga.exe',
#                     shell=True,stdout=open(os.devnull,"w"))
#             ll_i = self.glue.readlp_raw('lineinfo.dat')
#             ll = self.glue.con_lptonp_raw(ll_i)

#             # copy orginal linelist into a working list
#             self.ll = ll.copy()

#             # copy line list into ll.dat file
#             shutil.copy('lineinfo.dat','/dev/shm/FAL/{0}/ll.dat'.format(self.ID))

#             # delete old fort.11 file
#             if os.path.isfile("fort.11"):
#                 os.unlink('fort.11')
#             if os.path.isfile('/dev/shm/FAL/{0}/fort.11'.format(self.ID)):
#                 os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))

#             # write out the ascii line list
#             self.glue.writelp('/dev/shm/FAL/{0}/ll.dat'.format(self.ID),'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
#             os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')

#             # reinitialize SYNTHE Python Wrapper
#             SYNTHE.reset()
#             self.SYNTHE = synthepkg.synthe(ID=ID,verbose=False,clobber=True,rline=False,LL=True,rlinedict={},arct_bool=self.arct_bool)
#             self.SYNTHE.synbeg(self.indict,verbose=False)

#             print("Pro: {0} --> Full set of lines took {1}s".format(self.IDraw,time.time()-starttime))

#         else:
#             print("Pro: {0} --> Running User Defined Line List File...".format(self.IDraw))
#             # No need to reset SYNTHE
#             # set original SYNTHE object into self
#             self.SYNTHE = SYNTHE

#             # read in original linelist
#             orgll_i = self.glue.readlp(startll) 
#             self.orgll = self.glue.con_lptonp(orgll_i)
        
#             # copy orginal linelist into a working list
#             self.ll = self.orgll.copy()

#             # copy original line list into ll.dat file
#             if startll != '/dev/shm/FAL/{0}/ll.dat'.format(self.ID):
#                 shutil.copy(startll,'/dev/shm/FAL/{0}/ll.dat'.format(self.ID))
     
#         # read in output line list and make it the one that we are going to work on
#         newll =  self.ll.copy()
#         # print(np.unique(newll['CODE']))

#         # put a fail safe to not go through the gamma switch step
#         if linelistcheck:
#             os.chdir(self.parentdir)
#             return

#         # pull gamma info for only lines in self.ll, if it isn't found just set to gamma_W
#         self.gammaswitch = {}
#         self.gammaswitch['switch'] = np.array(['W' if (float(SPECODE) < 100) else 'M' for SPECODE in newll['CODE']])
#         self.gammaswitch['GSWITCH_str'] = np.array(
#             ['{0:7.4f}+{1:7.2f}'.format(wl,float(code)) for wl,code in zip(newll['WL'],newll['CODE'])]
#             )

#         for ii,gsw_str in enumerate(self.gammaswitch['GSWITCH_str']):
#             matind = np.in1d(gammafile["GSWIT_str"],gsw_str,assume_unique=True)
#             if any(matind):
#                 self.gammaswitch['switch'][ii] = gammafile['GAMMAPAR'][matind].data[0]
#             else:
#                 pass

#         # determine runtime speed
#         if "speed" in self.indict.keys():
#             self.speed = self.indict['speed']
#         else:
#             self.speed = None

#         # cd back into parent directory
#         os.chdir(self.parentdir)
        
#     def __call__(self,p,lineind=None,reusell=False,writespec=False,verbose=False,timeit=False,fullgamma=False,tau_i=False,speed=None,trans=None):
        
#         os.chdir('{0}'.format(self.ID))
        
#         '''
#         Define free parameters, stick them into the line list
#         and then call SYNTHE to get spectrum
#         '''
#         # start if timer if needed
#         starttime = time.time()
 
#         if speed != None:
#             self.speed = speed
            
#         # have condition where p = []
#         if len(p) == 0:
#             pass
#         else:
#             # delete old fort.11 file
#             os.unlink('fort.11')
#             os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))

#             # define free parameters and stick them into line list
#             # FREE POSIBLE PARS: DWL, DGFLOG, and one or all of DGAMMAR, DGAMMAS, DGAMMAW
    
#             # the case where all lines are free (i.e., no lineind index array)
#             if lineind==None:
#                 del self.ll['DWL','DGFLOG','DGAMMAR','DGAMMAS','DGAMMAW']
#                 self.ll.add_columns(p.columns.values())
#             else:
#                 # else just replace the parameters for the lines that are indictated in lineind
#                 for pind,llind in enumerate(lineind):
#                     self.ll['DWL'][llind] = p['DWL'][pind]
#                     self.ll['DGFLOG'][llind] = p['DGFLOG'][pind]
#                     if fullgamma:
#                         self.ll['DGAMMAR'][llind] = p['DGAMMAR'][pind]
#                         self.ll['DGAMMAS'][llind] = p['DGAMMAS'][pind]
#                         self.ll['DGAMMAW'][llind] = p['DGAMMAW'][pind]
#                     else:
#                         if self.gammaswitch['switch'][pind] == 'W':
#                             gammapar = "DGAMMAW"
#                         elif self.gammaswitch['switch'][pind] == 'R':
#                             gammapar = "DGAMMAR"
#                         elif self.gammaswitch['switch'][pind] == 'S':
#                             gammapar = "DGAMMAS"
#                         else:
#                             raise IOError('COULD NOT FIGURE OUT GAMMA SWITCH!!!')
#                         self.ll[gammapar][llind] = p['DGAMMA'][pind]

#             # convert the table into the correct string format
#             lpfmttab = self.glue.con_nptolp(self.ll)

#             # write out the ascii line list
#             self.glue.writelp(lpfmttab,'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
#             os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')

#         if reusell:
#             self.SYNTHE.reset()

#         if timeit:
#             print("Pro: {1} --> Adjusted p array -- Step time: {0:7.5f} s".format(time.time()-starttime,self.IDraw))
#             lasttime = time.time()

#         # now run the synthe steps
#         # -- read line list --
#         if (verbose == True or verbose == 'readlines'):
#             verbose_i = True
#         else:
#             verbose_i = False
#         self.SYNTHE.readlines(rtype='punch500',rlinedict={},verbose=verbose_i,reuse_ll=reusell)
#         if timeit:
#             print("Pro: {1} --> Read in lines -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#             lasttime = time.time()

#         # -- run xnfpelsyn --
#         if (verbose == True or verbose == 'xnfpelsyn'):
#             verbose_i = True
#         else:
#             verbose_i = False
#         self.SYNTHE.xnfpelsyn(verbose=verbose_i)
#         if timeit:
#             print("Pro: {1} --> XNFPELSYN -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#             lasttime = time.time()

#         # -- do synthe --
#         if (verbose == True or verbose == 'syn'):
#             verbose_i = True
#         else:
#             verbose_i = False
#         self.SYNTHE.syn(verbose=verbose_i,speed=self.speed)
#         if timeit:
#             print("Pro: {1} --> SYNTHE -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#             lasttime = time.time()

#         # -- do spectrv --
#         if (verbose == True or verbose == 'spectrv'):
#             verbose_i = True
#         else:
#             verbose_i = False
#         self.SYNTHE.spectrv(verbose=verbose_i,tau=tau_i)
#         if timeit:
#             print("Pro: {1} --> SPECTRV -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#             lasttime = time.time()

#         # -- do rotate --
#         if (verbose == True or verbose == 'rotate'):
#             verbose_i = True
#         else:
#             verbose_i = False
#         self.SYNTHE.rotate(self.indict,verbose=verbose_i)
#         if timeit:
#             print("Pro: {1} --> ROTATE -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#             lasttime = time.time()

#         # -- check if the user wants broadening --
#         if self.indict['MACVEL'] == -1:
#             # clean out fort.21 (input)
#             if os.path.exists('fort.21'):
#                 os.unlink('fort.21')
#             if os.path.isfile('/dev/shm/FAL/{0}/fort.21'.format(self.ID)):
#                 os.remove('/dev/shm/FAL/{0}/fort.21'.format(self.ID))
            
#             # mv input file into fort.21
#             self.SYNTHE._makesym('/dev/shm/FAL/{0}/{1}'.format(self.ID,'ROT1'),'ROT1_mac_inst')

#             print("Pro: {1} --> No Broadening Applied -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
    
#         else:
#             # -- do broaden for macroturblence --
#             if (verbose == True or verbose == 'broad_mac'):
#                 verbose_i = True
#             else:
#                 verbose_i = False
#             self.SYNTHE.broaden('ROT1',self.indict,broadtype='MAC',verbose=verbose_i)
#             if timeit:
#                 print("Pro: {1} --> BROADEN MACRO -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#                 lasttime = time.time()

#             # -- do broaden for instrumental --
#             if (verbose == True or verbose == 'broad_inst'):
#                 verbose_i = True
#             else:
#                 verbose_i = False
#             self.SYNTHE.broaden('ROT1_mac',self.indict,broadtype='INSTRUMENT',WLreg=self.WLreg,verbose=verbose_i)
#             if timeit:
#                 print("Pro: {1} --> BROADEN INSTR -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#                 lasttime = time.time()

#         # read in binary output spectrum
#         # self.outspec = self.glue.readspecbin('ROT1_mac_inst',N=int(10000000))
#         runinfo = {}
#         with open('RUNINFO.dat','r') as runinfoff:
#             runinfoarr = runinfoff.readlines()
#             runinfo_i = runinfoarr[1].split()
#             runinfo['LENGTH'] = int(runinfo_i[0])
#             runinfo['RATIO'] = float(runinfo_i[1])
#             runinfo['WBEGIN'] = float(runinfo_i[2])
#             runinfo['DWLBEG'] = float(runinfo_i[3])
#             runinfo['WLLAST'] = float(runinfo_i[4])
#             runinfo['DWLLAST'] = float(runinfo_i[5])

#         self.outspec = self.glue.readspecbin('ROT1_mac_inst',N=runinfo['LENGTH'])

#         # self.outspec = self.glue.readspecbin('ROT1')
#         if timeit:
#             print("Pro: {1} --> Read in binary spectrum -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#             lasttime = time.time()

#         if trans != None:
#             self.outspec = np.multiply(self.outspec['QMU1'],trans['FLUX'])

#             if timeit:
#                 print("Pro: {1} --> Multiply Transmission Spectrum -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#                 lasttime = time.time()


#         try:
#             assert self.outspec
#             # write out ascii version of spectrum, line list, and header file
#             if writespec:
#                 # remove any old ascii output files
#                 outfiles = ['headinfo.dat','lineinfo.dat','specfile.dat']
#                 for ofile in outfiles:
#                     try:
#                         os.remove(ofile)
#                     except OSError:
#                         pass
#                 os.unlink('fort.1')
#                 shutil.copy('/dev/shm/FAL/{0}/ROT1_mac_inst'.format(self.ID),'fort.1')
#                 subprocess.check_call('/home1/02349/cargilpa/FAL/PYTHON/bin/syntoascanga.exe',
#                     shell=True,stdout=open(os.devnull,"w"))
#             if timeit:
#                 print("Pro: {1} --> Write ASCII files -- Step time: {0:7.5f} s".format(time.time()-lasttime,self.IDraw))
#                 lasttime = time.time()

#         except (AssertionError,AttributeError):
#             self.outspec = None

#         os.chdir(self.parentdir)

#         if timeit:
#             print("Pro: {1} --> FINISHED FALmod -- Full time: {0:7.5f} s".format(time.time()-starttime,self.IDraw))
        
#         return self.outspec

#     def __del__(self):
#         os.chdir(self.parentdir)


#   WLin,DWLin,GFLOGin,DGFLOGin,CODEin,Ein,XJin,LABELin,&
#   EPin,XJPin,LABELPin,GRin,DGAMMARin,GSin,DGAMMASin,GWin,DGAMMAWin,WAVENOin,&
#   REFin,NBLOin,NBUPin,ISO1in,X1in,ISO2in,X2in,OTHER1in,OTHER2in,ISOSHIFTin,&
#   NELIONin,residin

