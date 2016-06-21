from __future__ import print_function
import sys, time
import FALmcmc as FALmcmc
import time, shutil, os, multiprocessing
import numpy as np
# import psutil
from astropy.table import Table
import traceback

import warnings
warnings.simplefilter(action='ignore',category=FutureWarning)


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

	print('... Setting up Seg{0}'.format(IDin))

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
		outputfile=outputfile)
	try:
		if testrun:
			return MCMC
		else:
			print("Seg{0} - Working on {1}".format(IDin,multiprocessing.current_process().name))

			# 100 walkers, 500 steps

			# # build samplers
			# MCMC.buildsampler(nwalkers=100,threads=0)
			# # run MCMC
			# MCMC.run_MCMC(500,burnin=False,nburn=0)

			# build samplers
			MCMC.buildsampler(nwalkers=150,threads=0)
			# run MCMC
			MCMC.run_MCMC(500,burnin=False,nburn=0)
			
			return 	"RUN TOOK {0:10.5f}".format(time.time()-MCMC.starttime)

	except Exception as e:
		print('Caught Exception in worker thread (Seg{0} - Working on {1})'.format(IDin,multiprocessing.current_process().name))
		traceback.print_exc()
		print()
		raise e

def makeinlist(infilename):
	# timing info
	starttime = time.time()
	walltime = time.time() + float(os.environ.get('SLURM_TACC_RUNLIMIT_MINS'))*60.0

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
	regfile = Table.read(infilename,names=['ID','WLstart','WLend','LINWLstart','LINWLend','WLRAN','NUMLINES','ARCSCALE'],format='ascii')

	indictlist = []

	for ii,rf_i in enumerate(regfile):
		tempdict = ({'starttime':starttime,'walltime':walltime,'IDin':int(rf_i['ID']),
			'minWLin':float(rf_i['WLstart']),'maxWLin':float(rf_i['WLend']),
			'minlinWL':float(rf_i['LINWLstart']),'maxlinWL':float(rf_i['LINWLend']),
			'arcscale':float(rf_i['ARCSCALE']),
			'outputfile':'MCMC_{0}.dat'.format(rf_i['ID']),
			'RUNID':ii})

		indictlist.append(tempdict)

	return indictlist

if __name__ == '__main__':
	# run makeinlist
	indictlist = makeinlist(sys.argv[1])

	print("READ IN {0} NUMBER OF SEGMENTS".format(len(indictlist)))

	MPI = False

	##### multiprocessing stuff #######
	pool_size = multiprocessing.cpu_count()
	# pool_size = 2*multiprocessing.cpu_count()
	# pool_size = 8

	os.system("taskset -pc 0-%d %s" % (pool_size,os.getpid()))

	if MPI:
		# Initialize the MPI-based pool used for parallelization.
		pool = MPIPool(debug=False)

		if not pool.is_master():
			# Wait for instructions from the master process.
			pool.wait()
			sys.exit(0)

	else:

		pool = multiprocessing.Pool(processes=pool_size)



	# # ONE RUN TEST:
	# resultprt = runFAL({'starttime':starttime,'walltime':walltime,'IDin':3,
	# 	'minWLin':596.0,'maxWLin':598.0,
	# 	'initlines':'/work/02349/cargilpa/FAL/DATA/SLpars.dat',
	# 	'condst':seg3cond,
	# 	'outputfile':'MCMC_5960_5980.dat'})
	# print resultprt

	# set up CPU precentage reader
	# psutil.cpu_percent(interval=1, percpu=True)


	print("FOUND POOL SIZE OF {0} CPUS".format(pool_size))

	pool.map(runFAL, indictlist)
	pool.close() # no more tasks
	pool.join()  # wrap up current tasks

	# print('... CPU Percentage Usage:')
	# print(psutil.cpu_percent(percpu=True))