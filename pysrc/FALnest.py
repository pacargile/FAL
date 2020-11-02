import sys,time,os
from datetime import datetime
import dynesty
import h5py

import numpy as np
from astropy.table import Table,vstack,unique
from scipy.interpolate import UnivariateSpline
from scipy.stats import beta
from scipy.optimize import minimize
from scipy.optimize import SR1
from scipy import constants
speedoflight = constants.c / 1000.0

import FALlinesel
import FALmod
import FALGlue
from FALGlue import *

import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)

conroypath = os.environ.get('CSCRATCH')+'pac'
holypath   = os.environ.get('HOLYSCRATCH')+'conroy_lab/pacargile'
homepath   = os.environ.get('HOME')

datapath = holypath

def lnprob(pin,args,verbose=False):

    # read in arguments
    (solobswave,solobsflux,arcobswave,arcobsflux,transflux,bg_sol_flux,bg_arc_flux,fmdict,fmll,Tarr,minWL,maxWL,minLWL,maxLWL) = args
    # for printing to command line, generate an ID
    fmdictkeys = list(fmdict)
    ID = fmdictkeys[0] - 10000000

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

    # arcturus scaling and vel shift
    arcscale = pin[-4]
    arcvel = pin[-3]

    # transmission spectrum scaling and vel shift
    transcale = pin[-2]
    tranvel   = pin[-1]

    # check to make sure transcale is > 0.0
    if transcale < 0.0:
        print('transcale issue {0}'.format(transcale))
        return -np.inf,[np.nan,np.nan]
    # apply gaussian prior to transscale
    tsprior = -0.5*( ((transcale-1.0)**2.0)/(0.001**2.0))

    if (tranvel < -1.0) or (tranvel > 1.0):
        print('tranvel issue {0}'.format(tranvel))
        return -np.inf,[np.nan,np.nan]

    # check arcturus scalings
    if (arcscale > 1.05) or (arcscale < 0.95):
        print('arcscale issue {}'.format(arcscale))
        return -np.inf,[np.nan,np.nan]
    if (arcvel < -1.0) or (arcvel > 1.0):
        print('arcvel issue {}'.format(arcvel))
        return -np.inf,[np.nan,np.nan]

    # scale transmission spectrum, renormalize it
    transflux_i = ((transcale*transflux)-(transcale-1.0))
    # shfit transmission spectrum and resample it back to obswave
    transflux_ii = UnivariateSpline(obswave['Sun']*(1.0+(tranvel/speedoflight)),transflux_i,s=0,k=1)(obswave['Sun'])

    # now divide the observed spectrum by the transmission spectrum
    obsflux['Sun'] = obsflux['Sun'] / transflux_ii

    # scale arcturus spectrum
    obsflux['Arcturus'] = obsflux['Arcturus']*arcscale
    obswave['Arcturus'] = obswave['Arcturus']*(1.0+(arcvel/speedoflight))

    # divide observed spectrum by background models
    obsflux['Sun'] = obsflux['Sun'] / bg_sol_flux
    obsflux['Arcturus'] = obsflux['Arcturus'] / bg_arc_flux

    # calculate lnl
    lnl, mod = lnlike(p,obswave,obsflux,fmdict,minWL,maxWL)

    # print(pin,gplp,lp,lnl)

    # modblob = mod

    # return the lnprob value
    lnpr = lnl
    return lnpr
    
def lnlike(p,obswave,obsflux,fmdict,minWL,maxWL):
    # generate ID list
    IDlist = list(fmdict)
    IDlist.sort()

    sig = {}

    sig['Sun'] = np.ones_like(obsflux['Sun'])/200.0
    sig['Arcturus'] = np.ones_like(obsflux['Arcturus'])/200.0

    # sig['Sun'] = obsflux['Sun']/500.0
    # sig['Arcturus'] = obsflux['Arcturus']/300.0

    # sig['Sun'] = np.ones_like(obsflux['Sun'])/1000.0
    # sig['Arcturus'] = np.ones_like(obsflux['Arcturus'])/1.0

    modintrp = {}

    # initialize lnp
    lnp = 0
    # loop over stars and calclate lnp_i
    for ID_i,star_i in zip(IDlist,['Sun','Arcturus']):
        # calculate model spectrum for p
        _spec,_ll = fmdict[ID_i].runsynthe(timeit=False,linelist='readlast',parr=p,verbose=False)
        _spectab_i = Table(_spec)
        _spectab = _spectab_i[(_spectab_i['WAVE'] <= maxWL) & (_spectab_i['WAVE'] >= minWL)]
        _specflux = _spectab['QMU1']/_spectab['QMU2']
        _specintr = UnivariateSpline(_spectab['WAVE'].data,_specflux,s=0,k=1,ext=1)(obswave[star_i])
        modintrp[star_i] = _specintr
        # residsq = (np.subtract(obsflux[star_i],_specintr)**2.0)/(sig[star_i]**2.0)
        # lnp_i = np.sum(-0.5*residsq + np.log(1.0/np.sqrt(2*np.pi*(sig[star_i]**2.0))))
        # lnp = lnp+lnp_i

        residsq = np.divide((np.subtract(obsflux[star_i],_specintr)**2.0),(sig[star_i]**2.0))
        residsqnorm = np.add(-0.5*residsq,np.log(1.0/np.sqrt(2*np.pi*(sig[star_i]**2.0))))
        lnp_i = np.sum(residsqnorm) 
        lnp = lnp+lnp_i

    return lnp, modintrp

def priortrans(upars,args):
    pflag = args

    pars = []
    for upars_i,pflag_i in zip(upars[:-4],pflag):
        if pflag_i == 'WL':
            pars_i = 0.02*upars_i - -0.01

        if pflag_i == 'GF':
            pars_i = (1.5 - -10.0)*upars_i - -10.0

        if pflag_i == 'GW':
            pars_i = (0.65 - -1.5)*upars_i - -1.5

        pars.append(pars_i)

    pars.append( (1.05 - 0.95) *upars[-4] + 0.95)
    pars.append( (1.0  - -1.0) *upars[-3] + -1.0)
    pars.append( (1.25 -  0.0) *upars[-2] + 0.0)
    pars.append( (1.0  - -1.0) *upars[-1] + -1.0)

    print(pars)

    return pars

class FALnest(object):
    def __init__(self,**kwargs):

        minWLin     = kwargs.get("minWLin",605.2)
        maxWLin     = kwargs.get("maxWLin",605.8)
        minlinWL    = kwargs.get("minlinWL",None)
        maxlinWL    = kwargs.get("maxlinWL",None)
        IDin        = kwargs.get("IDin",1)
        starttime   = kwargs.get("starttime",None)
        walltime    = kwargs.get("walltime",None)
        initlines   = kwargs.get("initlines",None)
        injectlines = kwargs.get("injectlines",None)
        outputfile  = kwargs.get("outputfile",None)
        outputdir  = kwargs.get("outputdir",None)
        arcscale    = kwargs.get("arcscale",None)

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
        self.wavebuffer = 0.15
        self.waverange = [self.minWL-self.wavebuffer,self.maxWL+self.wavebuffer]

        # setting cut for line selection
        self.condst = [{'LP':'RESID','OP':np.less,'LV':0.99}]

        # set output file name
        self.outputfile = outputfile
        self.outputdir = outputdir

        if arcscale == None:
            self.arcscale = 1.0
        else:
            self.arcscale = arcscale

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
        print('Seg: {0} --> Total possible run time = {1} seconds'.format(self.ID,self.walltime-self.starttime))
        print('Seg: {0} --> Full Wavelength Range = {1:9.5f}-{2:9.5f} ({3:9.5f}) nm'.format(
            self.ID,self.minWL-self.wavebuffer,self.maxWL+self.wavebuffer,self.maxWL-self.minWL+(2.0*self.wavebuffer)))


        # set up some dictionaries for passing objects
        self.fmdict = {}
        origsyndict = {}

        # run synthe to get all needed lines
        for ID_i,star_i in zip(self.IDlist,['Sun','Arcturus']):
            print('Seg: {0} --> Running original synthe on full master line list for {1}'.format(self.ID,star_i))
            # initialize the class
            fm_i = FALmod.FALmod(ID=ID_i,waverange=self.waverange,starpars=star_i,verbose=False)
            # run SYNTHE using the master line list to grab all important lines
            spec_i,ll_i = fm_i.runsynthe(timeit=False,linelist='readmaster')
            self.fmdict[ID_i] = fm_i
            origsyndict[ID_i] = [Table(spec_i),ll_i]

        # Assemble working line list (union of ll_i from last for-loop)
        print('Seg: {0} --> Assemble working line list'.format(self.ID))
        # stack the tables
        fmll = vstack([origsyndict[ID_i][1] for ID_i in self.IDlist])
        # fmll['FILTERBOOL'] = np.zeros(len(fmll['WL']),dtype=int)
        # sort tables on all collumns, include RESID as that way the stronger line will be listed first and will be set as a fit parameter
        tabpars = ([
            'WL','GFLOG', 'CODE', 'E', 'XJ', 'LABEL', 'EP', 'XJP', 'LABELP', 
            'GR', 'GS', 'GW', 'WAVENO', 'REF', 'NBLO', 'NBUP', 'ISO1', 'X1', 'ISO2', 'X2', 
            'OTHER1', 'OTHER2'])
        fmll.sort(tabpars+['RESID'])
        fmll = unique(fmll,tabpars)

        # check to see if working in H-band where we don't incldue a background model
        if ((self.minWL > 1300.0) & (self.maxWL < 2300.0)):
            print("Seg: {0} --> Not using a background model for weak lines".format(self.ID))
        else:
            # remove lines with RESID > 0.9999 as these are included in the static background model
            fmll = fmll[fmll['RESID'] < 0.9999]
            print("Seg: {0} --> Using a background model for weak lines".format(self.ID))

        # set it into self
        self.fmll = fmll

        # remove all predicted lines with depths < 0.99
        plinesind = np.array((self.fmll['LABEL'] == '          ') & (self.fmll['LABELP'] == '          ') & (self.fmll['RESID'] < 0.99),dtype=bool)
        self.fmll = self.fmll[~plinesind]

        # change all DWL back to zero, hack for pervious solar-only fit
        self.fmll['DWL'] = np.zeros_like(self.fmll['DWL'])

        # now run each stellar spectrum again and archive the results
        print('Seg: {0} --> Archiving the results into working directories'.format(self.ID))
        for ID_i,star_i in zip(self.IDlist,['Sun','Arcturus']):
            if star_i == 'Arcturus':
                print('Seg: {0} --> Changing C12/C13 abundance fraction for Arcturus'.format(self.ID))
                C12C13 = 7.4
                fmll_i = self.fmll.copy()
                for ll_ii in fmll_i:
                    if ll_ii['ISO1'] == 12:
                        ll_ii['GFLOG'] = ll_ii['GFLOG']+0.005-np.log10(1.0+(1.0/C12C13))
                    elif ll_ii['ISO1'] == 13:
                        ll_ii['GFLOG'] = ll_ii['GFLOG']+1.955-np.log10(1.0+C12C13)
                    else:
                        pass
            else:
                fmll_i = self.fmll.copy()
            fmll_i.remove_column('RESID')
            _spec,_ll = self.fmdict[ID_i].runsynthe(timeit=False,linelist=fmll_i,archive=True)
            # reset orgll to fmll because we don't want to use the synthe parsed ll
            self.fmdict[ID_i].orgll = fmll_i

        # run function to select which lines are modeled
        (self.parr,self.psig,self.pflag,self.Tarr) = FALlinesel.linesel(self.fmll,self.condst,self.minLINWL,self.maxLINWL)

        # calculate number of lines with free paramters and their wavelengths
        fitlineind = []
        fitlinewl  = []
        for ii,tt in enumerate(self.Tarr):
            if not all(tt == -1):
                fitlineind.append(ii)
                fitlinewl.append(float(self.fmll['WL'][ii]))


        print("Seg: {0} --> Total number of lines considered in WL segment = {1}".format(self.ID,len(self.fmll)))
        print("Seg: {0} --> Number of lines that are free in WL segment = {1}".format(self.ID,len(fitlinewl)))
        # print("Seg: {0} --> Index in line list of the modeled lines...".format(self.ID),fitlineind)
        # print("Seg: {0} --> WL of modeled lines...".format(self.ID),fitlinewl)

        # number of dimensions
        self.ndim = len(self.parr)
        print("Seg: {0} --> Number of Free Line Parameters...".format(self.ID),self.ndim)
        print("Seg: {0} --> Fitting Transmission Spectrum (scaling and velocity)".format(self.ID))
        self.ndim = self.ndim + 2
        print("Seg: {0} --> Fitting Arcturus Spectrum scaling and velocity".format(self.ID))
        self.ndim = self.ndim + 2

        # get observed data, transmission spectrum, and background model
        (self.solobswave,self.solobsflux,self.arcobswave,self.arcobsflux,self.transflux,self.bg_sol_flux,self.bg_arc_flux) = self.getspecdata()

        # # initialize output files
        # self.initoutput()

        # # now that the output files have been written...
        # # for arcturus, remove troublesome pixels in observed spectrum (flux < 0.001 & flux > 0.99)
        # prunearc = (self.arcobsflux >= 0.001) & (self.arcobsflux <= 0.99)
        # self.arcobsflux = self.arcobsflux[np.array(prunearc)]
        # self.arcobswave = self.arcobswave[np.array(prunearc)]

        self.arcobsflux = self.arcobsflux*self.arcscale

        print("Seg: {0} --> Number of Pixels in Obs Sol ...".format(self.ID),len(self.solobswave))
        print("Seg: {0} --> Number of Pixels in Obs Arc ...".format(self.ID),len(self.arcobswave))
        print("Seg: {0} --> Using an initial scaling of {1} for Obs Arc ...".format(self.ID,self.arcscale))

        print("Seg: {0} --> Finished Setup".format(self.ID))


    def getspecdata(self):
        if ((self.minWL > 470.0) & (self.maxWL < 755.0)):
            print("Seg: {0} --> Working with Blue Optical Spectrum".format(self.ID))
            # read in observed data for the Sun and Acturus
            sol_i = Table.read('{0}/FAL/data/SOL_4750_8270.fits'.format(datapath),format='fits')            
            arc_ii = Table.read('{0}/FAL/data/ARC_3800_9300_HINKLE.fits'.format(datapath),format='fits')
            arc_i = Table()
            arc_i['WAVE'] = arc_ii['WAVELENGTH']/10.0
            arc_i['FLUX'] = arc_ii['ARCTURUS'].copy()
            fluxcond = arc_i['FLUX'] >= 0.0
            fluxcond = np.array(fluxcond,dtype=bool)
            print("Seg: {0} --> Number of pixels clipped out of Arcturus spectrum: {1}".format(self.ID,len(arc_i)-len(fluxcond.nonzero()[0])))
            arc_i = arc_i[fluxcond]

            # read in transmission spectrum
            transh5 = h5py.File('{0}/FAL/data/TRANS_OPT_10_22_15.h5'.format(datapath),'r')
            trans = Table(np.array(transh5['spec']))
            trans.sort('WAVE')
            # correct for slight doppler shift
            trans['WAVE'] = trans['WAVE']*(1.0-(1.231/speedoflight))

            # read in background spectrum for each model
            bg_sol_i = Table.read('{0}/FAL/data/SPEC_SOL_weak_475_1000.fits.gz'.format(datapath),format='fits')
            bg_arc_i = Table.read('{0}/FAL/data/SPEC_ARC_weak_475_1000.fits.gz'.format(datapath),format='fits')


        elif ((self.minWL > 745.0) & (self.maxWL < 1010.0)):
            print("Seg: {0} --> Working with Red Optical Spectrum".format(self.ID))
            # read in observed data for the Sun and Acturus
            sol_i = Table.read('{0}/FAL/data/SOL_710_1098.fits'.format(datapath),format='fits')            
            arc_ii = Table.read('{0}/FAL/data/ARC_3800_9300_HINKLE.fits'.format(datapath),format='fits')
            arc_i = Table()
            arc_i['WAVE'] = arc_ii['WAVELENGTH']/10.0
            arc_i['FLUX'] = arc_ii['ARCTURUS'].copy()
            fluxcond = arc_i['FLUX'] >= 0.0
            fluxcond = np.array(fluxcond,dtype=bool)
            print("Seg: {0} --> Number of pixels clipped out of Arcturus spectrum: {1}".format(self.ID,len(arc_i)-len(fluxcond.nonzero()[0])))
            arc_i = arc_i[fluxcond]

            # read in transmission spectrum
            transh5 = h5py.File('{0}/FAL/data/TRANS_RED_2_21_17.h5'.format(datapath),'r')
            trans = Table(np.array(transh5['spec']))
            trans.sort('WAVE')
            # correct for slight doppler shift
            trans['WAVE'] = trans['WAVE']*(1.0-(1.231/speedoflight))

            # read in background spectrum for each model
            bg_sol_i = Table.read('{0}/FAL/data/SPEC_SOL_weak_475_1000.fits.gz'.format(datapath),format='fits')
            bg_arc_i = Table.read('{0}/FAL/data/SPEC_ARC_weak_475_1000.fits.gz'.format(datapath),format='fits')

        elif ((self.minWL > 1300.0) & (self.maxWL < 2300.0)):
            print("Seg: {0} --> Working with H-Band Spectrum".format(self.ID))
            sol_i  = Table.read('{0}/FAL/data/SOL_HBAND_Kur_8_26_15.fits'.format(datapath),format='fits')
            arc_ii = Table.read('{0}/FAL/data/ARC_HBAND_HINKLE.fits'.format(datapath),format='fits')
            arc_i = Table()
            arc_i['WAVE'] = (arc_ii['Wavelength_air'].copy()/10.0)*(1.0+(-12.1/speedoflight))
            arc_i['FLUX'] = arc_ii['Flux'].copy()
            fluxcond = arc_i['FLUX'] >= 0.0
            fluxcond = np.array(fluxcond,dtype=bool)
            print("Seg: {0} --> Number of pixels clipped out of Arcturus spectrum: {1}".format(self.ID,len(arc_i)-len(fluxcond.nonzero()[0])))
            arc_i = arc_i[fluxcond]

            # read in transmission spectrum
            transh5 = h5py.File('{0}/FAL/data/TRANS_HBAND_10_22_15.h5'.format(datapath),'r')
            trans = Table(np.array(transh5['spec']))
            trans.sort('WAVE')

            # no background model for H-band spectrum
            bg_sol_i = np.ones_like(sol_i['FLUX'])
            bg_arc_i = np.ones_like(arc_i['FLUX'])

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

        # check to see if we are working in the H-band with no background model spectrum
        if bg_sol_i != np.ones_like(sol_i['FLUX']):
            # parse and interpolate background model spectrum
            bg_sol = bg_sol_i[ (bg_sol_i['WAVE'] > solobswave.min()-0.1) & (bg_sol_i['WAVE'] < solobswave.max()+0.1) ]
            bg_sol['FLUX'] = bg_sol['QMU1']/bg_sol['QMU2']
            bg_sol_flux = UnivariateSpline(bg_sol['WAVE'].data,bg_sol['FLUX'].data,s=0,k=1)(solobswave)

            # check to see if there is observed ARC spectra in wavelength range
            if len(arcobswave) > 0:
                bg_arc = bg_arc_i[ (bg_arc_i['WAVE'] > arcobswave.min()-0.1) & (bg_arc_i['WAVE'] < arcobswave.max()+0.1) ]
                bg_arc['FLUX'] = bg_arc['QMU1']/bg_arc['QMU2']
                bg_arc_flux = UnivariateSpline(bg_arc['WAVE'].data,bg_arc['FLUX'].data,s=0,k=1)(arcobswave)
            else:
                bg_arc_flux = np.array([])

        else:
            bg_sol_flux = bg_sol_i
            bg_arc_flux = bg_arc_i

        return (solobswave,solobsflux,arcobswave,arcobsflux,transflux,bg_sol_flux,bg_arc_flux)

    def initoutput(self):
        print("Seg: {0} --> Initializing Output Files".format(self.ID))
        # set up outfile
        if self.outputfile == None:
            self.outputfile = 'MCMC_{0:n}.dat'.format(time.time())
        if self.outputdir == None:
            self.outputdir = './'

        with open(self.outputdir+self.outputfile, "w") as f:
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
            f.write('WN\t')
            for ii in range(self.ndim):
                f.write('par{0}\t'.format(ii))
            f.write('lnprob\n')

        # WRITE LL WITH PARR INFO
        if os.path.exists(self.outputdir+'LL'+self.outputfile[4:-3]+'fits'):
            os.remove(self.outputdir+'LL'+self.outputfile[4:-3]+'fits')
        self.ll_i = self.fmll.copy()
        self.ll_i['FWL']     = self.Tarr[:,0]
        self.ll_i['FGFLOG']  = self.Tarr[:,1]
        self.ll_i['FGAMMAW'] = self.Tarr[:,2]
        self.ll_i['FGAMMAR'] = self.Tarr[:,3]
        self.ll_i['FGAMMAS'] = self.Tarr[:,4]
        self.ll_i.write(self.outputdir+'LL'+self.outputfile[4:-3]+'fits',format='fits',overwrite=True) 

        # Initialize HDF5 File for SAMP & write Sun observed spec
        if os.path.exists(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5'):
            os.remove(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5')
        self.soloutspec = Table()
        self.soloutspec['SOL_WAVE'] = self.solobswave
        self.soloutspec['SOL_FLUX'] = self.solobsflux
        self.soloutspec.write(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='solobs',overwrite=True)

        # write arcturus obs spec
        self.arcoutspec = Table()
        self.arcoutspec['ARC_WAVE'] = self.arcobswave
        self.arcoutspec['ARC_FLUX'] = self.arcobsflux
        self.arcoutspec.write(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='arcobs',overwrite=True,append=True)

        # compute zero spectrum with all the previous shifts applied
        ogspecdict = {}
        for ID_i,star_i in zip(self.IDlist,['Sun','Arcturus']):
            _spec,_ll = self.fmdict[ID_i].runsynthe(timeit=False,linelist='readlast',parr=self.fmll['DWL','DGFLOG','DGAMMAR','DGAMMAS','DGAMMAW'])
            _spectab = Table(_spec)
            _spectab.write(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path=star_i+'_zero',overwrite=True,append=True)

        # write original transmission spectrum into file
        transpec = Table()
        transpec['FLUX'] = self.transflux
        transpec.write(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='trans',overwrite=True,append=True)

        # write background model into file
        backgmod_s = Table()
        backgmod_s['SOL_FLUX'] = self.bg_sol_flux
        backgmod_a = Table()
        backgmod_a['ARC_FLUX'] = self.bg_arc_flux
        backgmod_s.write(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='SOL_background',overwrite=True,append=True)
        backgmod_a.write(self.outputdir+'SAMP'+self.outputfile[4:-3]+'h5',format='hdf5',path='ARC_background',overwrite=True,append=True)

        return 

    def run_nest(self):
        res = self._nest()
        return res

    def _nest(self):
        # flush STDOUT before chain starts, just to pring info to log
        sys.stdout.flush()

        inargs = ([(self.solobswave,self.solobsflux,
            self.arcobswave,self.arcobsflux,
            self.transflux,
            self.bg_sol_flux,self.bg_arc_flux,
            self.fmdict,self.fmll,self.Tarr,
            self.minWL,self.maxWL,
            self.minLINWL,self.maxLINWL
            )])

        dysampler = dynesty.NestedSampler(
               lnprob,
               priortrans,
               self.ndim,
               logl_args=inargs,
               ptform_args=[self.pflag],
               nlive=60,
               bound='multi',
               sample='rwalk',
               bootstrap=0,
               walks=5,
               )

        dysampler.run_nested(dlogz=0.5)
        return dysampler