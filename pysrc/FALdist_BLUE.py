from __future__ import print_function
import sys, time
import FALmcmc as FALmcmc
import time, shutil, os, multiprocessing, argparse
import numpy as np
# import psutil
from astropy.table import Table
import traceback

import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)

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

def runFAL(indict):
    starttime = indict['starttime']
    walltime = indict['walltime']
    IDin = indict['IDin']
    minWLin = indict['minWLin']
    maxWLin = indict['maxWLin']
    minlinWL=indict['minlinWL']
    maxlinWL=indict['maxlinWL']
    arcscale = indict['arcscale']
    runID = indict['RUNID']

    print('Seg: {0} --> Setting up segment'.format(IDin))

    # test keys in indict
    if 'testrun' in indict.keys():
        if indict['testrun'] == True:
            testrun = True
        else:
            testrun = False
    else:
        testrun = False

    if 'outputfile' in indict.keys():
        outputfile = indict['outputfile']
    else:
        outputfile = None

    if 'outputdir' in indict.keys():
        outputdir = indict['outputdir']
    else:
        outputdir = None

    if 'initlines' in indict.keys():
        initlines = indict['initlines']
    else:
        initlines = None


    # remove old working directory if it exisits
    try:
        shutil.rmtree('./{0}'.format(str(IDin).rjust(8,str(0))))
    except OSError:
        pass
    try:
        shutil.rmtree('/dev/shm/FAL/{0}'.format(str(IDin).rjust(8,str(0))))
    except OSError:
        pass

    # sleep if running more than 8 jobs at one time, to keep from 
    # maxing out memory
    if runID > 8:
        time.sleep(300)

    # setup FALmcmc
    MCMC = FALmcmc.FALmcmc(
        minWLin=minWLin,maxWLin=maxWLin,
        minlinWL=minlinWL,maxlinWL=maxlinWL,
        IDin=IDin,
        starttime=starttime,walltime=walltime,
        arcscale=arcscale,
        outputfile=outputfile,
        outputdir=outputdir,
        initlines=initlines)
    try:
        if testrun:
            return MCMC
        else:
            print("Seg: {0} --> Working on {1}".format(IDin,multiprocessing.current_process().name))

            # 100 walkers, 500 steps

            # build samplers
            MCMC.buildsampler(nwalkers=100,threads=0)
            # run MCMC
            MCMC.run_MCMC(500,burnin=False,nburn=0)

            # # build samplers
            # MCMC.buildsampler(nwalkers=250,threads=0)
            # # run MCMC
            # MCMC.run_MCMC(3000,burnin=False,nburn=0)
            
            print("Seg: {0} --> RUN TOOK {1:10.5f}".format(IDin,time.time()-MCMC.starttime))
            return 

    except Exception as e:
        print('Caught Exception in worker thread (Seg: {0}, Working on {1})'.format(IDin,multiprocessing.current_process().name))
        traceback.print_exc()
        print()
        raise e

def makeinlist(infilename):
    # timing info
    starttime = time.time()
    try:
        walltime = time.time() + float(os.environ.get('SLURM_TACC_RUNLIMIT_MINS'))*60.0
    except TypeError:
        walltime = time.time() + 10000000000.0

    # set up input dictionary

    # list of conditional statements for line selection, the list may contain
    # multiple conditional statements for each segement
    #
    # format of conditional statements (must have at:
    # {'LP': ... , -> Name of line parameter
    #  'OP': ... , -> conditional operator, use numpy conditional functions
    #  'LV': ... } -> value you want to select LP on
    #
    # for numpy conditional statements, look at:
    # http://docs.scipy.org/doc/numpy/reference/routines.logic.html
    
    # read slicer file
    try:
        regfile = Table.read(infilename,names=['ID','WLstart','WLend','LINWLstart','LINWLend','WLRAN','NUMLINES'],format='ascii')
        regfile['ARCSCALE'] = np.ones_like(regfile['WLstart'])
    except:
        regfile = Table.read(infilename,names=['ID','WLstart','WLend','LINWLstart','LINWLend','WLRAN','NUMLINES','ARCSCALE'],format='ascii')

    indictlist = []

    MASTERLL = (['{0}/FAL/MASTERLL/FULLOPT/KuruczLL_450_1350.bin'.format(datapath),
        '{0}/FAL/MASTERLL/TiO/TiO_LL_450_1300.bin'.format(datapath)])

    for ii,rf_i in enumerate(regfile):
        tempdict = ({'starttime':starttime,'walltime':walltime,'IDin':int(rf_i['ID']),
            'minWLin':float(rf_i['WLstart']),'maxWLin':float(rf_i['WLend']),
            'minlinWL':float(rf_i['LINWLstart']),'maxlinWL':float(rf_i['LINWLend']),
            'arcscale':float(rf_i['ARCSCALE']),
            'outputfile':'MCMC_{0}.dat'.format(rf_i['ID']),
            'outputdir':'/n/holyscratch01/conroy_lab/pacargile/FAL/test/output/',
            'initlines':None,
            # 'initlines':'/n/regal/conroy_lab/pac/FAL/OPTSEG/CarRun1_LL_475_750_new.h5',
            # 'initlines':MASTERLL,
            'RUNID':ii})

        indictlist.append(tempdict)

    return indictlist

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--segfile',
        help='Name of file containing segment information',
        type=str,default='./SEG_OPT.dat')
    parser.add_argument('-i','--index',
        help='Index in segment file for which one user wants to run',
        type=int,default=0)

    args = parser.parse_args()

    # run makeinlist
    indictlist = makeinlist(args.segfile)

    if args.index <= len(indictlist):
        indict = indictlist[args.index]
        print("Seg: {0} --> Read in segment with index {1} from {2}".format(
            indict['IDin'],args.index,args.segfile))
        runFAL(indict)
    else:
        print("NOT ENOUGH ROWS IN SEG FILE FOR THIS RUN!")
