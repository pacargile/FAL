from __future__ import print_function
import sys,time
from datetime import datetime
import emcee
import h5py

import numpy as np
from astropy.table import Table,vstack,unique
from scipy.interpolate import UnivariateSpline
from scipy.stats import beta
from scipy import constants
speedoflight = constants.c / 1000.0

import FALlinesel
import FALmod
import FALGlue
from FALGlue import *

def lnprob(pin,args,verbose=False,justprior=False):

    # read in arguments
    (solobswave,solobsflux,arcobswave,arcobsflux,transflux,fmdict,fmll,Tarr,minWL,maxWL) = args
    # for printing to command line, generate an ID
    ID = fmdict.keys()[0] - 10000000

    # set up observed data dictionaries
    obswave = {'Sun':solobswave,'Arcturus':arcobswave}
    obsflux = {'Sun':solobsflux,'Arcturus':arcobsflux}

    # stick new p into Table
    p = Table()
    p['DWL']     = np.array([pin[xx] if xx != -1 else fmll['DWL'][kk]     for kk,xx in enumerate(Tarr[...,0])])
    p['DGFLOG']  = np.array([pin[xx] if xx != -1 else fmll['DGFLOG'][kk]  for kk,xx in enumerate(Tarr[...,1])])
    p['DGAMMAW'] = np.array([pin[xx] if xx != -1 else fmll['DGAMMAW'][kk] for kk,xx in enumerate(Tarr[...,2])])
    p['DGAMMAR'] = np.array([pin[xx] if xx != -1 else fmll['DGAMMAR'][kk] for kk,xx in enumerate(Tarr[...,3])])
    p['DGAMMAS'] = np.array([pin[xx] if xx != -1 else fmll['DGAMMAR'][kk] for kk,xx in enumerate(Tarr[...,4])])

    # transmission spectrum scaling and vel shift
    transcale = pin[-2]
    tranvel   = pin[-1]

    # check to make sure transcale is > 0.0
    if transcale < 0.0:
        return -np.inf,[np.nan,np.nan]
    if (tranvel < -1.0) or (tranvel > 1.0):
        return -np.inf,[np.nan,np.nan]

    # calcluate prior
    lp = lnprior(p,ID,Tarr,fmll,minWL,maxWL,verbose=verbose)
    try:
        if np.isfinite(lp) == False:
                return -np.inf, [np.nan,np.nan]
    except ValueError:
        pass
    # give out for just calculating the prior probability
    if justprior:
        return lp

    # scale transmission spectrum, renormalize it
    transflux_i = ((transcale*transflux)-(transcale-1.0))
    # shfit transmission spectrum and resample it back to obswave
    transflux_ii = UnivariateSpline(obswave['Sun']*(1.0+(tranvel/speedoflight)),transflux_i,s=0,k=1)(obswave['Sun'])

    # now divide the observed spectrum by the transmission spectrum
    obsflux['Sun'] = obsflux['Sun'] / transflux_ii

    # calculate lnl
    lnl, mod = lnlike(p,obswave,obsflux,fmdict,minWL,maxWL)

    # print(pin,gplp,lp,lnl)

    modblob = mod

    # return the lnprob value
    lnpr = lp + lnl
    return lnpr, modblob
    
def lnlike(p,obswave,obsflux,fmdict,minWL,maxWL):
    # generate ID list
    IDlist = fmdict.keys()
    IDlist.sort()

    # sigma based on SNR = 1000
    sig = {}
    sig = {'Sun':1.0/1000.0,'Arcturus':1.0/300.0}

    modintrp = {}

    # initialize lnp
    lnp = 0
    # loop over stars and calclate lnp_i
    for ID_i,star_i in zip(IDlist,['Sun','Arcturus']):
        # calculate model spectrum for p
        _spec,_ll = fmdict[ID_i].runsynthe(timeit=False,linelist='readlast',parr=p)
        _spectab_i = Table(_spec)
        _spectab = _spectab_i[(_spectab_i['WAVE'] <= maxWL) & (_spectab_i['WAVE'] >= minWL)]
        _specflux = _spectab['QMU1']/_spectab['QMU2']
        _specintr = UnivariateSpline(_spectab['WAVE'].data,_specflux,s=0,k=1,ext=1)(obswave[star_i])
        modintrp[star_i] = _specintr
        residsq = (np.subtract(obsflux[star_i],_specintr)**2.0)/(sig[star_i]**2.0)
        lnp_i = np.sum(-0.5*residsq + np.log(1.0/np.sqrt(2*np.pi*(sig[star_i]**2.0))))
        lnp = lnp+lnp_i

    return lnp, modintrp



def lnprior(p,ID,Tarr,fmll,minWL,maxWL,verbose=False):

    # apply non-informative prior on wavelength to make sure
    # line is not shifted outside working segment
    for ii,pp in enumerate(zip(p['DWL'][Tarr[...,0] != -1],fmll['WL'][Tarr[...,0] != -1])):
        if (np.abs(pp[0]) > 0.0) & (pp[1] > minWL) & (pp[1] < maxWL):
            wlshift = pp[1]+pp[0]
            if (wlshift < minWL-0.05) or (wlshift > maxWL+0.05):
                if verbose:
                    print('Pro: {0} --> CAUGHT A WAVELENGTH SHIFT OUTSIDE SPECTRUM BOUNDS, LINE SHIFTED TO: {1}'.format(ID,wlshift))
                return -np.inf

    # Prior on gamma using beta function
    mingamma = -1.0
    maxgamma = 1.0
    rangegamma = maxgamma-mingamma
    gammawprior = beta.logpdf(p['DGAMMAW'][Tarr[...,2] != -1],1.0,1.0,loc=mingamma,scale=rangegamma)
    gammarprior = beta.logpdf(p['DGAMMAR'][Tarr[...,3] != -1],1.0,1.0,loc=mingamma,scale=rangegamma)
    gammasprior = beta.logpdf(p['DGAMMAS'][Tarr[...,4] != -1],1.0,1.0,loc=mingamma,scale=rangegamma)

    # check to see if it returns any priors outside uniform prior
    if (any(np.isinf(gammawprior)) or 
        any(np.isinf(gammarprior)) or 
        any(np.isinf(gammasprior))) :
        if verbose:
            print('Pro: {0} --> CAUGHT A GAMMA SHIFT OUTSIDE THE PRIORS'.format(ID))
        return -np.inf

    # Prior on gf using beta function
    mingflog = -10.0
    maxgflog = 2.0
    rangegflog = maxgflog-mingflog
    gfprior = beta.logpdf(p['DGFLOG'][Tarr[...,1] != -1],1.0,1.0,loc=mingflog,scale=rangegflog)

    # check to see if it returns any priors outside uniform prior
    if any(np.isinf(gfprior)):
        if verbose:
            print('Pro: {0} --> CAUGHT A LOG(GF) SHIFT OUTSIDE THE PRIORS'.format(ID))
        return -np.inf

    velshift = 75.0 #km/s
    wsh = fmll['WL'][Tarr[...,0] != -1]*(velshift/speedoflight)
    wsh_max = max(wsh)
    minwll = -1.0*wsh_max
    maxwll = wsh_max
    rangewll = maxwll-minwll
    wlprior = beta.logpdf(p['DWL'][Tarr[...,0] != -1],2.0,2.0,loc=minwll,scale=rangewll)

    # check to see if it returns any priors outside uniform prior
    if any(np.isinf(wlprior)):
        if verbose:
            for ii,wlp in enumerate(wlprior):
                if np.isinf(wlp):
                    print('Pro: {0} --> CAUGHT A WAVELENGTH SHIFT OUTSIDE THE PRIORS: line {2} shifted by {1} nm'.format(ID,p['DWL'][ii],fmll['WL'][ii]))
        return -np.inf

    # check to see if arrays are empty, if so add 0.0 so that the sumation works
    if len(gammawprior) == 0:
        gammawprior = [0.0]
    if len(gammarprior) == 0:
        gammarprior = [0.0]
    if len(gammasprior) == 0:
        gammasprior = [0.0]

    # 2-D gaussian prior on delta(log(gf)) and delta(lambda)
    # reparameterize such that they move on space evenly
    sig_dgflog = 0.5
    sig_dWL = 1.0
    gf_wl_prior = (-0.5*((p['DGFLOG'][Tarr[...,1] != -1]/sig_dgflog)**2.0)*((p['DWL'][Tarr[...,0] != -1]/sig_dWL)**2.0)) #- 0.5*np.log(2.0*np.pi*(sig_coup**2.0))

    # RETURN WITH COUPLED PRIOR
    return np.sum(np.hstack([gammawprior,gammarprior,gammasprior,gfprior,wlprior,gf_wl_prior]))

class FALmcmc(object):
    def __init__(self,minWLin=605.2,maxWLin=605.8,
        minlinWL=None,maxlinWL=None,
        IDin=1,
        starttime=None,walltime=None,
        initlines=None,injectlines=None,
        outputfile=None):

        # change the following lines to be inputs when you init the class
        self.minWL = minWLin
        self.maxWL = maxWLin
        if minlinWL == None:
            self.minLINWL = self.minWL
            self.maxLINWL = self.maxWL
        else:                
            self.minLINWL = minlinWL
            self.maxLINWL = maxlinWL

        self.ID = IDin
        self.numstars = 2
        self.IDlist = [int(10000000*x)+self.ID for x in range(1,self.numstars+1,1)]

        # Define synthesis wavelength range
        self.waverange = [self.minWL-0.1,self.maxWL+0.1]

        # setting cut for line selection
        self.condst = [{'LP':'RESID','OP':np.less,'LV':0.99}]

        # set output file name
        self.outputfile = outputfile

        # define a general start time so that the code can stop short of the wall time and save everything
        if starttime == None:
            self.starttime = time.time()
        else:
            self.starttime = starttime # starttime in seconds

        if walltime == None:
            # an impossibly large walltime
            self.walltime = self.starttime + 300000
        else:
            self.walltime = walltime # walltime in seconds
        print('Pro: {0} --> Total possible run time = {1} seconds'.format(self.ID,self.walltime-self.starttime))
        print('Pro: {0} --> Full Wavelength Range = {1:9.5f}-{2:9.5f} ({3:9.5f}) nm'.format(self.ID,self.minWL-0.1,self.maxWL+0.1,self.maxWL-self.minWL+0.2))


        # set up some dictionaries for passing objects
        self.fmdict = {}
        origsyndict = {}

        # run synthe to get all needed lines
        for ID_i,star_i in zip(self.IDlist,['Sun','Arcturus']):
            print('Pro: {0} --> Running original synthe on full master line list for {1}'.format(self.ID,star_i))
            # initialize the class
            fm_i = FALmod.FALmod(ID=ID_i,waverange=self.waverange,starpars=star_i)
            # run SYNTHE using the master line list to grab all important lines
            spec_i,ll_i = fm_i.runsynthe(timeit=False,linelist='readmaster')
            self.fmdict[ID_i] = fm_i
            origsyndict[ID_i] = [Table(spec_i),ll_i]

        # Assemble working line list (union of ll_i from last for-loop)
        print('Pro: {0} --> Assemble working line list'.format(self.ID))
        # stack the tables
        fmll = vstack([origsyndict[ID_i][1] for ID_i in self.IDlist])
        # fmll['FILTERBOOL'] = np.zeros(len(fmll['WL']),dtype=int)
        # sort tables on all collumns, include RESID as that way the stronger line will be listed first and will be set as a fit parameter
        tabpars = ['WL','GFLOG', 'CODE', 'E', 'XJ', 'LABEL', 'EP', 'XJP', 'LABELP', 'GR', 'GS', 'GW', 'WAVENO', 'REF', 'NBLO', 'NBUP', 'ISO1', 'X1', 'ISO2', 'X2', 'OTHER']
        fmll.sort(tabpars+['RESID'])
        fmll = unique(fmll,tabpars)

        # # inject fake lines
        # fmll = self.injectfake(fmll.copy())

        # # inject previous parameters
        # fmll = self.injectprev(fmll.copy())

        # set it into self
        self.fmll = fmll

        # now run each stellar spectrum again and archive the results
        print('Pro: {0} --> Archiving the results into working directories'.format(self.ID))
        for ID_i in self.IDlist:
            _spec,_ll = self.fmdict[ID_i].runsynthe(timeit=False,linelist=self.fmll,archive=True)
            # reset orgll to fmll because we don't want to use the synthe parsed ll
            self.fmdict[ID_i].orgll = self.fmll

        # run function to select which lines are modeled
        (self.parr,self.psig,self.pflag,self.Tarr) = FALlinesel.linesel(self.fmll,self.condst,self.minLINWL,self.maxLINWL)

        # calculate number of lines with free paramters and their wavelengths
        fitlineind = []
        fitlinewl  = []
        for ii,tt in enumerate(self.Tarr):
            if not all(tt == -1):
                fitlineind.append(ii)
                fitlinewl.append(float(self.fmll['WL'][ii]))


        print("Pro: {0} --> Total number of lines considered in WL segment = {1}".format(self.ID,len(self.fmll)))
        print("Pro: {0} --> Number of lines that are free in WL segment = {1}".format(self.ID,len(fitlinewl)))
        # print("Pro: {0} --> Index in line list of the modeled lines...".format(self.ID),fitlineind)
        # print("Pro: {0} --> WL of modeled lines...".format(self.ID),fitlinewl)

        # number of dimensions
        self.ndim = len(self.parr)
        print("Pro: {0} --> Number of Free Line Parameters...".format(self.ID),self.ndim)
        print("Pro: {0} --> Fitting Transmission Spectrum (scaling and velocity)".format(self.ID))
        self.ndim = self.ndim + 2

        # get observed data and transmission spectrum
        (self.solobswave,self.solobsflux,self.arcobswave,self.arcobsflux,self.transflux) = self.getspecdata()

        # initialize output files
        self.initoutput()

        print("Pro: {0} --> Finished Setup".format(self.ID))

    def injectfake(self):
        # # print out if the code found any lines to inject
        # FAKELINES = 0
        # for ll_t in self.fmll:
        #     if 'FAK' in str(ll_t['REF']):
        #         print("Pro: {0} --> Inject FAKE line at WL={1}".format(self.ID,ll_t['WL']))
        #         FAKELINES = FAKELINES + 1
        # if FAKELINES == 0:
        #     print("Pro: {0} --> NO FAKE LINES FOUND IN SEGMENT".format(self.ID))
        return

    def injectprev(self):
        # # determine if there are any pre-initialized lines (from previous run) and set those free parameters
        # if initlines != None:
        #     # make a unique ID for each line
        #     fmllcode = np.empty(len(self.fm.ll),dtype=object)
        #     for nnn,ill in enumerate(self.fm.ll):
        #         fmllcode[nnn] = "".join(
        #             [str(ill['WL']),
        #             ill['CODE'],
        #             ill['E'],ill['EP'],
        #             ill['XJ'],ill['XJP'],
        #             ill['REF'],
        #             ill['ISO1'],ill['X1'],
        #             ill['ISO2'],ill['X2'],
        #             ill['OTHER']]
        #             ).replace(" ","")

        #     # Read previous table: LINE INFO, DWL, DGFLOG, DGAMMA (will figure out which GAMMA after the fact)
        #     ilines = Table(np.array(h5py.File(initlines,'r')['data']))

        #     ilines.sort('WL')
        #     # ilines = Table.read(initlines,format='fits')
        #     ilines['UNIQ_ID'] = np.empty(len(ilines),dtype=object)
        #     for nnn,ill in enumerate(ilines):
        #         ilines['UNIQ_ID'][nnn] = "".join(
        #             [str(ill['WL']),
        #             ill['CODE'],
        #             ill['E'],ill['EP'],
        #             ill['XJ'],ill['XJP'],
        #             ill['REF'],
        #             ill['ISO1'],ill['X1'],
        #             ill['ISO2'],ill['X2'],
        #             ill['OTHER']]
        #         ).replace(" ","")
        #     numpreset = 0
        #     for ii,fmlc in enumerate(fmllcode):
        #         cond_intl = np.in1d(ilines['UNIQ_ID'],fmlc,assume_unique=True)
        #         if any(cond_intl):
        #             numpreset = numpreset + 1
        #             # print("Pro: {0} --> Setting Previous Pars for WL = {1:7.4f}".format(self.ID,float(self.fm.ll['WL'][ii])))
        #             self.fm.ll['DWL'][ii] = float('{0:6.4f}'.format(float(ilines['DWL'][cond_intl])))
        #             self.fm.ll['DGFLOG'][ii] = float('{0:6.4f}'.format(float(ilines['DGFLOG'][cond_intl])))
        #             if self.fm.gammaswitch['switch'][ii] == 'W':
        #                 self.fm.ll['DGAMMAW'][ii] = float('{0:6.4f}'.format(float(ilines['DGAMMA'][cond_intl])))
        #             elif self.fm.gammaswitch['switch'][ii] == 'R':
        #                 self.fm.ll['DGAMMAR'][ii] = float('{0:6.4f}'.format(float(ilines['DGAMMA'][cond_intl])))
        #             elif self.fm.gammaswitch['switch'][ii] == 'S':
        #                 self.fm.ll['DGAMMAS'][ii] = float('{0:6.4f}'.format(float(ilines['DGAMMA'][cond_intl])))
        #             else:
        #                 pass
        #     print("Pro: {0} --> Setting Previous Pars for {1} lines".format(self.ID,numpreset))

        return

    def getspecdata(self):
        if ((self.minWL > 470.0) & (self.maxWL < 800.0)):
            print("Pro: {0} --> Working with Optical Spectrum".format(self.ID))
            # read in observed data for the Sun and Acturus
            sol_i = Table.read('/work/02349/cargilpa/FAL/DATA/SOL_4750_8270.fits',format='fits')            
            arc_ii = Table.read('/work/02349/cargilpa/FAL/DATA/ARC_3800_9300_HINKLE.fits',format='fits')
            arc_i = Table()
            arc_i['WAVE'] = arc_ii['WAVELENGTH']
            arc_i['FLUX'] = arc_ii['ARCTURUS']

            # read in transmission spectrum
            transh5 = h5py.File('/work/02349/cargilpa/FAL/DATA/TRANS/TRANS_OPT_10_22_15.h5','r')
            trans = Table(np.array(transh5['spec']))
            trans.sort('WAVE')
            # correct for slight doppler shift
            trans['WAVE'] = trans['WAVE']*(1.0-(1.231/speedoflight))

        elif ((self.minWL > 1300.0) & (self.maxWL < 2300.0)):
            print("Pro: {0} --> Working with H-Band Spectrum".format(self.ID))
            sol_i = Table.read('/work/02349/cargilpa/FAL/DATA/SOL_HBAND_Kur_8_26_15.fits',format='fits')
            arc_ii = Table.read('/work/02349/cargilpa/FAL/DATA/ARC_HBAND_HINKLE.fits',format='fits')
            arc_i = Table()
            arc_i['WAVE'] = arc_ii['Wavelength_air']
            arc_i['FLUX'] = arc_ii['Flux']
            arc_i.sort('WAVE')

            # read in transmission spectrum
            transh5 = h5py.File('/work/02349/cargilpa/FAL/DATA/TRANS/TRANS_HBAND_10_22_15.h5','r')
            trans = Table(np.array(transh5['spec']))
            trans.sort('WAVE')

        else:
            raise IOError("COULD NOT FIND OBSERVED SOLAR SPECTRUM IN THIS WAVELENGTH RANGE!")

        # parse solar spectrum
        SW = sol_i['WAVE'].data
        SF = sol_i['FLUX'].data
        solind = np.argwhere( (SW > self.minWL) & (SW < self.maxWL) )
        solobswave = np.squeeze(np.array(SW[solind]))
        solobsflux = np.squeeze(np.array(SF[solind]))

        # parse arcturus spectrum
        AW = arc_i['WAVE'].data
        AF = arc_i['FLUX'].data
        arcind = np.argwhere( (AW > self.minWL) & (AW < self.maxWL) )
        arcobswave = np.squeeze(np.array(AW[arcind]))
        arcobsflux = np.squeeze(np.array(AF[arcind]))

        # parse and interpolate transmission spectrum
        trans_i = trans[ (trans['WAVE'] > solobswave.min()-0.1) & (trans['WAVE'] < solobswave.max()+0.1) ]
        if 'FLUX' not in trans_i.keys():
            trans_i['FLUX'] = trans_i['QMU1'] / trans_i['QMU2']
        transintrp = UnivariateSpline(trans_i['WAVE'].data,trans_i['FLUX'].data,s=0,k=1)(solobswave)
        transflux = transintrp

        return (solobswave,solobsflux,arcobswave,arcobsflux,transflux)

    def initoutput(self):
        print("Pro: {0} --> Initializing Output Files".format(self.ID))
        # set up outfile
        if self.outputfile == None:
            self.outputfile = 'MCMC_{0:n}.dat'.format(time.time())

        with open(self.outputfile, "w") as f:
            f.write(
                '# RUN DATE/TIME: {0}, NumFP: {1} ,'.format(
                    datetime.today(),self.ndim))
            try:
                assert self.condst
                for condict in self.condst:
                    f.write('{0}'.format(condict['LP']))
                    if (condict['OP'].__str__() == "<ufunc 'less'>"):
                        f.write('<')
                    elif (condict['OP'].__str__() == "<ufunc 'greater'>"):
                        f.write('>')
                    elif (condict['OP'].__str__() == "<ufunc 'less_equal'>"):
                        f.write('<=')
                    elif (condict['OP'].__str__() == "<ufunc 'greater_equal'>"):
                        f.write('>=')
                    elif (condict['OP'].__str__() == "<ufunc 'equal'>"):
                        f.write('==')
                    elif (condict['OP'].__str__() == "<ufunc 'not_equal'>"):
                        f.write('!=')
                    else:
                        f.write('***')
                    f.write('{0} '.format(condict['LV']))
            except NameError:
                pass
            f.write('\n')
            f.write('# WAVELENGTH RANGE: {0:8.4f} {1:8.4f}'.format(self.minWL,self.maxWL))
            f.write(' LINE RANGE: {0:8.4f} {1:8.4f}\n'.format(self.minLINWL,self.maxLINWL))
            f.write('# TYPE OF FREE PARS: ')
            for pf in self.pflag:
                f.write('{0} '.format(pf))
            f.write('\n')
            f.write('WN ')
            for ii in range(self.ndim):
                f.write('par{0}\t'.format(ii))
            f.write('lnprob\n')

        # WRITE LL WITH PARR INFO
        self.ll_i = self.fmll.copy()
        self.ll_i['FWL']     = self.Tarr[:,0]
        self.ll_i['FGFLOG']  = self.Tarr[:,1]
        self.ll_i['FGAMMAW'] = self.Tarr[:,2]
        self.ll_i['FGAMMAR'] = self.Tarr[:,3]
        self.ll_i['FGAMMAS'] = self.Tarr[:,4]
        self.ll_i.write('LL'+self.outputfile[4:-3]+'fits',format='fits',overwrite=True) 

        # Initialize HDF5 File for SAMP & write Sun observed spec
        self.soloutspec = Table()
        self.soloutspec['SOL_WAVE'] = self.solobswave
        self.soloutspec['SOL_FLUX'] = self.solobsflux
        self.soloutspec.write('SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='solobs',overwrite=True)

        # write arcturus obs spec
        self.arcoutspec = Table()
        self.arcoutspec['ARC_WAVE'] = self.arcobswave
        self.arcoutspec['ARC_FLUX'] = self.arcobsflux
        self.arcoutspec.write('SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='arcobs',overwrite=True,append=True)

        # compute zero spectrum with all the previous shifts applied
        ogspecdict = {}
        for ID_i,star_i in zip(self.IDlist,['Sun','Arcturus']):
            _spec,_ll = self.fmdict[ID_i].runsynthe(timeit=False,linelist='readlast',parr=self.fmll['DWL','DGFLOG','DGAMMAR','DGAMMAS','DGAMMAW'])
            _spectab = Table(_spec)
            _spectab.write('SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path=star_i+'_zero',overwrite=True,append=True)

        # write original transmission spectrum into file
        transpec = Table()
        transpec['FLUX'] = self.transflux
        transpec.write('SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='trans',overwrite=True,append=True)

        return 

    def buildsampler(self,nwalkers=0,threads=0):
        if nwalkers == 0:
            # calculate number of walkers (default: 5 x number of free parameters)
            self.nwalkers = 5 * self.ndim
        else:	
            # use user defined number of walkers
            self.nwalkers = nwalkers

        print('Pro: {0} --> Number of walkers being used: {1}'.format(self.ID,self.nwalkers))

        # build sampler object
        args = ([(self.solobswave,self.solobsflux,
            self.arcobswave,self.arcobsflux,
            self.transflux,
            self.fmdict,self.fmll,self.Tarr,
            self.minWL,self.maxWL)])
        self.sampler = emcee.EnsembleSampler(
            self.nwalkers,self.ndim,
            lnprob,
            args=args,
            threads=threads,
            live_dangerously=True)

        # get p0 array and check for all finite values
        print('Pro: {0} --> Get inital walker positions'.format(self.ID))
        while True:
                # self.p0 = emcee.utils.sample_ball(self.parr,self.psig,self.nwalkers)
                self.p0 = self.buildball()
                testlp = [lnprob(pp,args[0],verbose=True,justprior=True) for pp in self.p0]
                if any(np.isinf(testlp)):
                    print('Pro: {0} --> ---- Need to redo p0 calculation'.format(self.ID))
                    print('Pro: {0} --> Problematic PAR: '.format(self.ID))
                    for li, lp_i in enumerate(testlp):
                        if np.isinf(lp_i):
                            print(self.p0[li],lp_i)
                    print('Pro: {0} --> ---- Reducing initial ball size by a factor of 1/2'.format(self.ID))
                    self.psig = self.psig * 0.5                        
                else:
                    print('Pro: {0} --> Initial Ball has following ranges...'.format(self.ID))
                    ballmin = np.amin(self.p0,axis=0)
                    ballmax = np.amax(self.p0,axis=0)
                    for ii in range(self.ndim):
                        print('Pro: {0} --> Par {1}: min = {2}, max = {3}'.format(self.ID,ii,ballmin[ii],ballmax[ii]))
                    break

    def buildball(self):
        p0out = []
        scalefact = 1.0
        velshift = 0.01 #km/s

        for _ in range(self.nwalkers):
            temparr = []            
            for ii,p in enumerate(self.parr):
                ps = self.psig[ii]
                pf = self.pflag[ii]

                if pf == 'WL':
                    fmll_i = self.ll_i[self.ll_i['FWL'] == ii]
                    if len(fmll_i) > 1:
                        wsh = max( fmll_i['WL'].data*(velshift/speedoflight) )
                        minwll = -1.0*wsh
                        maxwll = wsh
                        if (min( (fmll_i['WL']+fmll_i['DWL'])+minwll) < self.minWL) or (max( (fmll_i['WL']+fmll_i['DWL'])+maxwll) > self.maxWL):
                            minoff = min(
                                [min([np.abs(x+minwll-self.minWL) for x in (fmll_i['WL']+fmll_i['DWL'])]),
                                min([np.abs(x+maxwll-self.maxWL) for x in (fmll_i['WL']+fmll_i['DWL'])])]
                                )
                            minwll = -1.0*minoff
                            rangewll = 2.0*minoff
                        else:
                            rangewll = maxwll-minwll

                    else:
                        wsh = (fmll_i['WL']*(velshift/speedoflight))
                        minwll = -1.0*wsh
                        maxwll = wsh
                        if ( (fmll_i['WL']+fmll_i['DWL'])+minwll < self.minWL) or ((fmll_i['WL']+fmll_i['DWL'])+maxwll > self.maxWL):
                            minoff = min([np.abs( (fmll_i['WL']+fmll_i['DWL'])+minwll-self.minWL),np.abs( (fmll_i['WL']+fmll_i['DWL'])+maxwll-self.maxWL)])
                            minwll = -1.0*minoff
                            rangewll = 2.0*minoff
                        else:
                            rangewll = maxwll-minwll
                    wlshift = float(beta.rvs(2.0,2.0,loc=(minwll+fmll_i['DWL'][0])*scalefact,scale=rangewll*scalefact))
                    if (wlshift+fmll_i['WL'][0] < self.minWL-0.05) or (wlshift+fmll_i['WL'][0] > self.maxWL+0.05):
                        wlshift = np.zeros_like(fmll_i['DWL'][0]) + 0.0001*np.random.randn()
                    temparr.append(wlshift)

                elif pf == "GF":
                    fmll_i = self.ll_i[self.ll_i['FGFLOG'] == ii]
                    mingflog = -0.1
                    maxgflog = 0.1
                    rangegflog = maxgflog-mingflog
                    gflogshift = beta.rvs(4.0,4.0,loc=(mingflog+fmll_i['DGFLOG'][0])*scalefact,scale=rangegflog*scalefact)
                    if gflogshift >= 1.75:
                        # gflogshift = fmll_i['DGFLOG'][0] + 0.00001 * np.random.randn()
                        # gflogshift = np.zeros_like(fmll_i['DGFLOG'][0]) + 1.0*np.random.randn()
                        gflogshift = 1.25*np.ones_like(fmll_i['DGFLOG'][0]) + 0.1*np.random.randn()
                    if gflogshift <= -10.0:
                        gflogshift = -7.0*np.ones_like(fmll_i['DGFLOG'][0]) + 0.1*np.random.randn()

                    # temparr.append(beta.rvs(4.0,4.0,loc=(mingflog+fmll_i['DGFLOG'][0])*scalefact,scale=rangegflog*scalefact))
                    temparr.append(gflogshift)

                elif pf in ['GW','GS','GR']:
                    mingamma = -0.01
                    maxgamma = 0.01
                    rangegamma = maxgamma-mingamma
                    if pf == 'GW':
                        fmll_i = self.ll_i[self.ll_i['FGAMMAW'] == ii]
                        offsetgamma = fmll_i['DGAMMAW'][0]
                    elif pf == 'GS':
                        fmll_i = self.ll_i[self.ll_i['FGAMMAS'] == ii]
                        offsetgamma = fmll_i['DGAMMAS'][0]
                    elif pf == 'GR':
                        fmll_i = self.ll_i[self.ll_i['FGAMMAR'] == ii]
                        offsetgamma = fmll_i['DGAMMAR'][0]
                    else:
                        print('Pro: {0} --> ---- PROBELM WITH SETTING GAMMA OFFSET'.format(self.ID))
                        offsetgamma = 0.0
                    gammashift = beta.rvs(4.0,4.0,loc=(mingamma+offsetgamma)*scalefact,scale=rangegamma*scalefact)
                    if gammashift >= 0.9:
                        gammashift = 0.1*np.ones_like(offsetgamma) + 0.0001*np.random.randn()
                    if gammashift <= -0.9:
                        gammashift = -0.1*np.ones_like(offsetgamma) + 0.0001*np.random.randn()
                    temparr.append(gammashift)
                else:
                    pass

            # scaling for transmission spectrum
            temparr.append(0.1*np.random.randn()+1.0)
            # scaling for transmission velocity shift
            temparr.append(0.001*np.random.randn())

            # append the array to P0
            p0out.append(temparr)

        return p0out


    def run_MCMC(self,nsteps,burnin=True,nburn=100):
        if burnin:
            self.p0_i,self.lnprob0_i,self.rstate0_i = self._mcmc_bn(self.sampler,nburn,self.p0)
        else:
            self.p0_i = self.p0
            self.lnprob0_i = None
            self.rstate0_i = None
            
        self._mcmc(self.sampler,nsteps,self.p0_i,lnprob0=self.lnprob0_i,rstate0=self.rstate0_i)

    def _mcmc(self,sampler,niter,p0,lnprob0=None,rstate0=None):
        # RUN FINAL CHAIN MCMC
        print("Pro: {0} --> Starting MCMC with {1:n} links and {2:n} walkers".format(self.ID,niter,self.nwalkers))
        text = (
            "\rPro: {ID} --> MCMC step: {STEP:n}\n"+
            "Pro: {ID} --> Sampler's average acceptance fraction: {AF:5.2f}                  \n"+
            "Pro: {ID} --> Std(ln(Pr))/Avg(ln(Pr)) of Walker ln(Pr): {STAT:n}"
            )

        # get a walker number array
        walkernum = np.squeeze(np.arange(1,self.nwalkers+1,1).reshape(self.nwalkers,1))

        ii = 1
        # set print frequency
        prf = 1
        lasttime = time.time()
        deltaT = 0.0

        # define output file to append
        outf = open(self.outputfile,'a')

        outspec = h5py.File('SAMP'+self.outputfile[4:-3]+'h5','a')

        # set flag for final print
        finalflag = 0

        for pos,prob,_,blob in sampler.sample(
            p0,lnprob0=lnprob0,rstate0=rstate0,iterations=niter,storechain=False
            ):
            # set up output arrays
            pos_matrix=pos.reshape(self.nwalkers,self.ndim)
            # pos_matrix=np.append(walkernum,pos_matrix,axis=1)
            prob_array=np.squeeze(prob.reshape(self.nwalkers,1))
            steparray =np.column_stack([walkernum,pos_matrix,prob_array])
            # Write the current position to a file, one line per walker
            outf.write("\n".join(["\t".join([str(q) for q in p]) for p in steparray]))
            outf.write("\n")

            print(blob)

            outspec.create_dataset('SUN_{0}'.format(ii),data=blob['Sun'],compression='gzip')
            outspec.create_dataset('ARC_{0}'.format(ii),data=blob['Arcturus'],compression='gzip')

    		# handle SIGURS1 signal as a command to dump output file
            # try:
            #     assert GOT_SIG == False
            #     pass
            # except AssertionError:
            #     outf.flush()
            #     GOT_SIG = False

            # flush the buffer every X iterations
            if ((ii % 15 == 0.0) or (ii == niter)):
                outf.flush()
                sys.stdout.flush()
                outspec.flush()
            # # WRITE RESIDUAL TO OUTSPEC EVERY X iterations
            # if ((ii % 15 == 0.0) or (ii == niter)):
            #     outspec.flush()

            # print how the chains are doing
            #if ((ii % prf == 0.0) or (ii == niter)):
            print(
                text.format(
                    ID=self.ID,
                    STEP=ii,
                    AF=np.average(sampler.acceptance_fraction),
                    # STAT=(np.std(prob[goodind])/np.abs(np.average(prob[goodind])))
                    STAT=(np.nanstd(prob)/np.abs(np.nanmean(prob)))
                    )
                )
            print("Pro: {0} --> Step time: {1:5.3f}".format(self.ID,time.time()-lasttime))
            lasttime = time.time()
            ii = ii + 1
            # check if time is within ~1 min of wall time, if so stop the for-loop
            if ( (lasttime-self.starttime) > (self.walltime-120.0) ):
                print('Pro: {0} --> Stopping 2 min before walltime'.format(self.ID))
                finalflag = 1
		outf.flush()
		outspec.flush()
		outf.close()
		outspec.close()
                break
        if finalflag == 0:
            print('Pro: {0} --> Stopping at max iterations'.format(self.ID))            
            print('Pro: {0} --> Run time: {1} seconds'.format(self.ID,time.time()-self.starttime))
        outf.close()
        outspec.close()
