from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
from astropy.table import Table
import FALsynthepkg
import FALGlue
import shutil, os, subprocess
import numpy as np
import sys,glob,time,shutil
import broaden

import socket
hostname = socket.gethostname()
if hostname[:4] == 'holy':
     RUNLOC = 'ODY'
else:
     RUNLOC = 'LOCAL'

conroypath = os.environ.get('CSCRATCH')+'pac'
holypath   = os.environ.get('HOLYSCRATCH')+'conroy_lab/pacargile'
homepath   = os.environ.get('HOME')

datapath = holypath

import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)

class FALmod(object):
    """
    Class to take delta line parameters and return a synthesized
    spectrum from SYNTHE to use in python
    """    
    # def __init__(self,ID=None,starpars=None,waverange=None,verbose=False,extra={}):
    def __init__(self,**kwargs):
        '''
        FALmod -> Python wrapper around SYNTHE

        starpars dictionary can hold:
            - OBJECT: Sun or Arcturus
            - VROT: Rotation velocity (negative to turn on projection)
            - MACVEL: Macroturbulent velocity
            - OUTRES: Output resolution of spectrum (if OBJECT not Sun or Arcturs)
        '''


        # set ID
        ID = kwargs.get('ID',None)
        if ID is None:
            import uuid
            self.ID = uuid.uuid4().hex[:8]
            self.IDraw = None
        else:
            self.IDraw = ID
            self.ID = str(ID).rjust(8,str(0))

        # set starpars if starpars='Sun' or starpars='Arcturus'
        starpars = kwargs.get('starpars',None)
        self.starpars = {}
        if starpars is 'Sun':
            self.starpars['VROT'] = -2.02
            self.starpars['MACVEL'] = 1.5
            # self.starpars['MACVEL'] = -1
            self.starpars['OBJECT'] = 'Sun'
            self.starpars['RESOL'] = 3000000.0
        elif starpars is 'Arcturus':
            self.starpars['VROT'] = 1.0
            self.starpars['MACVEL'] = 3.0
            self.starpars['OBJECT'] = 'Arcturus'
            self.starpars['RESOL'] = 500000.0
        elif starpars is 'Mdwarf':
            self.starpars['VROT'] = 0.1
            self.starpars['MACVEL'] = -1.0
            self.starpars['OBJECT'] = 'Mdwarf'
            self.starpars['RESOL'] = 3000000.0
        else:
            if starpars is None:
                self.starpars = {'OBJECT':None}
            else:
                self.starpars = starpars

        waverange = kwargs.get('waverange',None)
        if waverange is None:
            self.starpars['WSTART'] = 1500.0
            self.starpars['WEND'] = 1501.0
        else:
            self.starpars['WSTART'] = waverange[0]
            self.starpars['WEND'] = waverange[1]

        extra = kwargs.get('extra',{})
        for kk in extra.keys():
            self.starpars[kk] = extra[kk]

        runinmem = kwargs.get('runinmem',True)
        if runinmem:
            # make dir in memory if not already there
            if not os.path.exists('/dev/shm/FAL/{0}'.format(self.ID)):
                os.makedirs('/dev/shm/FAL/{0}'.format(self.ID))

        verbose = kwargs.get('verbose',False)

        # make local dir to work in
        if not os.path.exists('{0}'.format(self.ID)):
            os.makedirs('{0}'.format(self.ID))

        # define parent directory
        self.parentdir = os.getcwd()

        # change into working directory
        os.chdir('{0}'.format(self.ID))

        # create the glue  
        self.glue = FALGlue.glue()

        # create broadening class
        self.brd = broaden.broaden()

        # initialize SYNTHE package
        self.SYNTHE = FALsynthepkg.synthe(ID=ID,verbose=verbose,clobber=True,starpars=self.starpars)

        # cd back into parent directory
        os.chdir(self.parentdir)


    def runsynthe(self,**kwargs):
        """
        Call SYNTHE
        """
        # sort out keywords
        # timing information printed out
        self.timeit = kwargs.get('timeit',False)
        # type of line list reader
        linelist = kwargs.get('linelist','readmaster')
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
        # print out the optical depth information
        self.tau_i = kwargs.get('tau_i',False)
        # definte master line list
        self.masterll = kwargs.get('masterll',None)

        # change into working directory
        os.chdir('{0}'.format(self.ID))

        # start if timer if needed
        if self.timeit:
            self.starttime = time.time()
            self.lasttime = self.starttime

        # call synbeg if readline != readlast

        if linelist is 'readlast':
            if not os.path.exists('/dev/shm/FAL/{0}/INT'.format(self.ID)):
                raise IOError("Pro: {1} --> WARNING: COULD NOT FIND INITIAL FILES!!!! {0:7.5f} s".format(time.time()-self.starttime,self.IDraw))
            else:
                self.SYNTHE.reset()
            # elif str(linelist) != 'readlast':
        else:
            if (verbose == True or verbose == 'synbeg'):
                verbose_i = True
            else:
                verbose_i = False
            self.SYNTHE.synbeg(self.starpars,clobber=True,verbose=verbose_i)
        # else:
        #     pass

        if archive:
            self.SYNTHE.archive()

        # adjust any delta values if parr has been provided
        if (verbose == True or verbose == 'adjustlines'):
            verbose_i = True
        else:
            verbose_i = False
        if parr != None:
            if str(linelist) == 'readlast':
                llo = self._adjustpar(parr,ll=self.orgll,verbose=verbose_i)
            else:
                llo = self._adjustpar(parr,verbose=verbose_i)

        # read in line list
        if (verbose == True or verbose == 'readlines'):
            verbose_i = True
        else:
            verbose_i = False
        self._readline(linelist,verbose_i)

        # do synthesis calc
        if type(verbose) == type(True):
            if (verbose == True):
                verbose_i = True
            else:
                verbose_i = False
        else:
            if ('synthesis' == verbose):
                verbose_i = verbose
            else:
                verbose_i = False
        self._synthesis(verbose=verbose_i)

        # do broadening
        if isinstance(verbose,bool):
            if (verbose == True):
                verbose_i = True
            else:
                verbose_i = False
        else:
            if ('broaden' == verbose):
                verbose_i = verbose
            else:
                verbose_i = False
        outspec,newll,binspecname = self._broaden(verbose=verbose)

        if archive:
            # write newll into INT/ for archiving purposes
            if os.path.exists('fort.11'):
                os.unlink('fort.11')
            if os.path.exists('/dev/shm/FAL/{0}/fort.11'.format(self.ID)):
                os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))
            newll_st = self.glue.con_nptolp(newll)
            self.glue.writelp(newll_st,'/dev/shm/FAL/{0}/INT/fort.11'.format(self.ID))
            self.orgll = newll.copy()

        if transspec != False:
            # do multiply transmission spectrum
            outspec = self._multiplytrans(outspec,transspec)

        if writespec:
            # write out spec, line list, and header to ascii
            self._writeout(binspecname)

        # cd back into parent directory
        os.chdir(self.parentdir)

        if self.timeit:
            print("Pro: {1} --> FINISHED FALmod -- Full time: {0:7.5f} s".format(time.time()-self.starttime,self.IDraw))

        return (outspec,newll)

    def __del__(self):
        os.chdir(self.parentdir)

    def _readline(self,linelist,verbose_i=False,speed=''):
        # read the line list appropriate for the user call
        # catch the case where the line list is a numpy table
        if type(linelist).__name__ == 'Table':
            # linelist equal to a path to user defined line list 
            self.SYNTHE.readlines(rtype=linelist,verbose=verbose_i)
            if speed == '':
                self.speed = 'fast'
            else:
                self.speed = speed

            if self.timeit:
                print("Pro: {1} --> Read in user defined numpy table line list -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()
            return

        elif linelist == 'readall':
            # read all individual line lists
            # rlinedict = {"TiO":True} # atoms, molecules + H2O & TiO
            # rlinedict = {"atoms":True,"moles":True,"H2O":True,"TiO":True,"predict":True} # atoms, molecules + H2O & TiO
            rlinedict = {"atoms":True,"moles":False,"H2O":False,"TiO":False,"predict":False} 
            self.SYNTHE.readlines(rtype='readall',rlinedict=rlinedict,verbose=verbose_i)
            if speed == '':
                self.speed = 'slow'
            else:
                self.speed = speed

            if self.timeit:
                print("Pro: {1} --> Read in all line lists -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()
            return

        elif linelist == 'readmaster':
            # determine which readmaster line list to use base on wavelength range
            if self.masterll == None:
                if (self.starpars['WSTART'] > 450.0) & (self.starpars['WEND'] < 1300.0):
                    MASTERLL = (['{0}/FAL/MASTERLL/FULLOPT/KuruczLL_450_1350.bin'.format(datapath),
                        '{0}/FAL/MASTERLL/TiO/TiO_LL_450_1300.bin'.format(datapath)])
                    # MASTERLL = (['/n/conroyfs1/pac/MASTERLL/FULLOPT/LL/Kur_LL_450_1300_CLEAN.bin'])
                elif (self.starpars['WSTART'] > 1399.0) & (self.starpars['WEND'] < 1901.0):
                    MASTERLL = (['{0}/FAL/MASTERLL/HBAND/KuruczLL_1400_1900.bin'.format(datapath),
                    '{0}/FAL/MASTERLL/H2O/KuruczH2OLL_1400_1900.bin'.format(datapath)])

                else:
                    print(self.starpars['WSTART'],self.starpars['WEND'])
                    raise ValueError('DID NOT UNDERSTAND MASTERLL')
            else:
                MASTERLL = self.masterll

            # read the masterline lists (cargile or kurucz plus H2O+TiO)
            self.SYNTHE.readlines(rtype='readmaster',verbose=verbose_i,MASTERLL=MASTERLL)
            if speed == '':
                self.speed = 'slow'
            else:
                self.speed = speed

            if self.timeit:
                print("Pro: {1} --> Read in master line list: {2} -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw,MASTERLL))
                self.lasttime = time.time()
            return

        elif linelist == 'readlast':
            # read the previously created line list in directory
            self.SYNTHE.readlines(rtype='readlast',verbose=verbose_i)
            if speed == '':
                self.speed = 'fast'
            else:
                self.speed = speed

            if self.timeit:
                print("Pro: {1} --> Read in lines -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()
            return

        else:
            # linelist equal to a path to user defined line list 
            self.SYNTHE.readlines(rtype=linelist,verbose=verbose_i)
            if speed == '':
                self.speed = 'fast'
            else:
                self.speed = speed

            linelistname = linelist
            if self.timeit:
                print("Pro: {1} --> Read in user defined line list {2} -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw,linelistname))
                self.lasttime = time.time()
            return


    def _adjustpar(self,parr_i,ll=None,verbose=False):
        # check to see if parr is one table of parameters (for the whole line list)
        # or if it is a dictionary of parameters and line indices
        if isinstance(parr_i,dict):
            parr = parr_i['parr']
            lineind = parr_i['lineind']
        else:
            parr = parr_i
            lineind = None

        # check to see if line list has been defined in self
        if ll == None:
            ll_i = self.glue.readlp_raw('/dev/shm/FAL/{0}/fort.11'.format(self.ID))
            ll = self.glue.con_lptonp_raw(ll_i)

        # copy orginal linelist into a working list
        llo = ll.copy()

        # delete old fort.11 file
        os.unlink('fort.11')
        os.remove('/dev/shm/FAL/{0}/fort.11'.format(self.ID))

        # define free parameters and stick them into line list
        # FREE POSIBLE PARS: DWL, DGFLOG, and one or all of DGAMMAR, DGAMMAS, DGAMMAW

        # the case where all lines are free (i.e., no lineind index array)
        if lineind==None:
            for kk in parr.keys():
                del llo[kk]
                llo[kk] = parr[kk]
        else:
            # else just replace the parameters for the lines that are indictated in lineind
            for pind,llind in enumerate(lineind):
                llo['DWL'][llind] = parr['DWL'][pind]
                llo['DGFLOG'][llind] = parr['DGFLOG'][pind]
                llo['DGAMMAR'][llind] = parr['DGAMMAR'][pind]
                llo['DGAMMAS'][llind] = parr['DGAMMAS'][pind]
                llo['DGAMMAW'][llind] = parr['DGAMMAW'][pind]

        # convert the table into the correct string format
        lpfmttab = self.glue.con_nptolp(llo)

        # write out the ascii line list
        self.glue.writelp(lpfmttab,'/dev/shm/FAL/{0}/fort.11'.format(self.ID))
        os.symlink('/dev/shm/FAL/{0}/fort.11'.format(self.ID),'fort.11')

        return llo

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
        self.SYNTHE.syn(verbose=verbose,speed=self.speed)
        if self.timeit:
            print("Pro: {1} --> SYNTHE -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

        # -- do spectrv --
        if (verbose == True or verbose == 'synthesis:spectrv'):
            verbose_i = True
        else:
            verbose_i = False
        self.SYNTHE.spectrv(verbose=verbose,tau=self.tau_i)
        if self.timeit:
            print("Pro: {1} --> SPECTRV -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

    def _broaden(self,verbose=False):
        # -- do rotate --
        if (verbose == True or verbose == 'broaden:rotate'):
            verbose_i = True
        else:
            verbose_i = False
        self.SYNTHE.rotate(self.starpars['VROT'],verbose=verbose)
        if self.timeit:
            print("Pro: {1} --> ROTATE -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()

        # pull outspec and newll from rotate.for code
        outspec,newll = self._specout('/dev/shm/FAL/{0}/{1}'.format(self.ID,'ROT1'),verbose=verbose)

        # -- do broaden for stellar broadening --
        if (verbose == True or verbose == 'broaden:broad_star'):
            verbose_i = True
        else:
            verbose_i = False

        # -- check if the user wants broadening --
        if 'MACVEL' in self.starpars.keys():
            if self.starpars['MACVEL'] == -1:
                # no macrovel applied, just return rotated spectrum to output
                print("Pro: {1} --> No Stellar Broadening Applied -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()
            else:
                vmacdict = {'type':'MACRO','units':'KM','val':self.starpars['MACVEL']}
                QMU1 = self.brd.broaden(outspec['WAVE'],outspec['QMU1'],vmacdict)
                outspec['QMU1'] = QMU1['FLUX']
                if self.timeit:
                    print("Pro: {1} --> BROADEN MACRO -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                    self.lasttime = time.time()

        # -- do broaden for instrumental --
        if (verbose == True or verbose == 'broaden:broad_inst'):
            verbose_i = True
        else:
            verbose_i = False

        if 'OUTRES' in self.starpars.keys():
            if self.starpars['OUTRES'] == -1:
                # no instrument broadening applied, just return rotated spectrum to output
                print("Pro: {1} --> No Instrument Broadening Applied -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()        
            else:
                instdict = {'type':'GAUSSIAN','units':'RESOLUTION','val':self.starpars['OUTRES']}
                QMU1 = self.brd.broaden(outspec['WAVE'],outspec['QMU1'],instdict)
                outspec['QMU1'] = QMU1['FLUX']

        else:
            if self.starpars['OBJECT'] == 'Sun': 
                if (self.starpars['WSTART'] > 450.0) & (self.starpars['WEND'] < 1300.0):               
                    instpars1 = self.SYNTHE.instparstr['OPT']['SINC']
                    instpars2 = self.SYNTHE.instparstr['OPT']['GAUSSIAN']
                elif (self.starpars['WSTART'] > 1399.0) & (self.starpars['WEND'] < 1901.0):
                    instpars1 = self.SYNTHE.instparstr['HBAND']['SINC']
                    instpars2 = self.SYNTHE.instparstr['HBAND']['GAUSSIAN']
                # the special case of the solar line profile with a SINC and a gaussian
                instdict1 = {'type':'SINX/X',  'val':instpars1,'units':'CM-1'}
                instdict2 = {'type':'GAUSSIAN','val':instpars2,'units':'CM-1'}

                for instdict_i in [instdict1,instdict2]:
                    QMU1 = self.brd.broaden(outspec['WAVE'],outspec['QMU1'],instdict_i)
                    outspec['QMU1'] = QMU1['FLUX']

            elif self.starpars['OBJECT'] == 'Arcturus':
                if (self.starpars['WSTART'] > 450.0) & (self.starpars['WEND'] < 1300.0):               
                    instpars = self.SYNTHE.instparstr['OPT']['GAUSSIAN']
                elif (self.starpars['WSTART'] > 1399.0) & (self.starpars['WEND'] < 1901.0):
                    instpars = self.SYNTHE.instparstr['HBAND']['GAUSSIAN']

                instdict = {'type':'GAUSSIAN','units':'RESOLUTION','val':instpars}
                QMU1 = self.brd.broaden(outspec['WAVE'],outspec['QMU1'],instdict)
                outspec['QMU1'] = QMU1['FLUX']
            elif self.starpars['OBJECT'] == 'Mdwarf':
                instdict = {'type':'GAUSSIAN','units':'RESOLUTION','val':instpars}
                QMU1 = self.brd.broaden(outspec['WAVE'],outspec['QMU1'],instdict)
                outspec['QMU1'] = QMU1['FLUX']
            else:
                # no instrument broadening applied, just return rotated spectrum to output
                print("Pro: {1} --> No Instrument Broadening Applied -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
                self.lasttime = time.time()

        if self.timeit:
            print("Pro: {1} --> BROADEN INSTR -- Step time: {0:7.5f} s".format(time.time()-self.lasttime,self.IDraw))
            self.lasttime = time.time()
        return (outspec,newll,'ROT1_mac_inst')

    def _specout(self,infile,verbose=False):
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
