import numpy as np
import os

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

def selmod(starpars):

	if starpars['OBJECT'] == 'Arcturus':
		# set up some object specific strings
		synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
			"0.          0   {LINOUT:3.0f}{TOL:7.5f}     {PRED}    00\n"
			"AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD")
		gausspar_opt = 130000.0
		# gausspar_hband = 100000.0
		gausspar_hband = 100000.0
		instparstr = {"OPT":{"GAUSSIAN":gausspar_opt},"HBAND":{"GAUSSIAN":gausspar_hband}}

		# modatm = os.path.expandvars("$HOME")+"/FAL/PYTHON/data/Arcturus_NEWpars_V4.dat"
		modatm = '{0}/FAL/ATM/arcturus_at12.atm'.format(datapath)

	elif starpars['OBJECT'] == 'Sun':
		# set synbeg string to use
		synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
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

		# set instrument profile pars
		instparstr = {"OPT":{},"HBAND":{}}
		sincpar_opt = 0.047628472
		gausspar_opt = 0.071442708
		sincpar_hband = 0.008140581
		gausspar_hband = 0.012210871
		instparstr['OPT'] = {"SINC":sincpar_opt,"GAUSSIAN":gausspar_opt}
		instparstr['HBAND'] = {"SINC":sincpar_hband,"GAUSSIAN":gausspar_hband}

		# set model atmosphere to use
		# modatm = os.path.expandvars("$HOME")+"/FAL/PYTHON/data/modcaspf.dat"
		# modatm = '/work/02349/cargilpa/FAL/TESTING/YSmod/YSsol.atm'
		modatm = '{0}/FAL/ATM/solar_at12.atm'.format(datapath)

	elif starpars['OBJECT'] == 'Mdwarf':
		# set up some object specific strings
		synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
			"0.          0    {LINOUT:3.0f} {TOL:7.5f}    {PRED}     0\n"
			"AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD")
		instparstr = {'OPT':None,'HBAND':None}

		instparstr = {'OPT':{"GAUSSIAN":150000.0},'HBAND':{'GAUSSIAN':150000.0}}
		# modatm = '{}/pac/FAL/data/ATM/mdwarf_at12.atm'.format(conroypath)
		modatm = '{0}/FAL/ATM/mdwarf_at12.atm'.format(datapath)

	else:
		# set up some object specific strings
		synbegvar = ("AIR       {WSTART:7.3f}   {WEND:7.3f}  {RESOL:10.1f} "
			"0.          0    {LINOUT:3.0f} {TOL:7.5f}    {PRED}     0\n"
			"AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD")
		instparstr = {'OPT':{"GAUSSIAN":100000.0},'HBAND':{'GAUSSIAN':100000.0}}

		# set model atmosphere to use **USING SOLAR BECAUSE THIS WILL CHANGE TO AN INTERPOLATOR**
		# modatm = os.path.expandvars("$HOME")+"/FAL/PYTHON/data/modcaspf.dat"
		# modatm = '{}/pac/FAL/data/ATM/arcturus_at12.atm'.format(conroypath)
		modatm = '{0}/FAL/ATM/arcturus_at12.atm'.format(datapath)

	return (synbegvar,instparstr,modatm)
