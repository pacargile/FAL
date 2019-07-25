from ctypes import cdll, CDLL, POINTER, ARRAY, c_int, c_double, c_char_p, c_char, c_void_p, c_float, c_long
import numpy as np
import fortranformat
from astropy.io import ascii
from astropy.table import Table
import subprocess
import copy

class glue(object):
	def __init__(self):

		try:
			# read in fortran libraray
			self.fortran = cdll.LoadLibrary('/n/conroyfs1/pac/FAL/lib/FALglue.so')
		except OSError:
			self.fortran = cdll.LoadLibrary('/Users/pcargile/Astro/FAL/PYSCRIPT/lib/FALglue.so')


		# define some useful things for later
		self.c_double_p = POINTER(c_double)
		self.c_float_p = POINTER(c_float)
		self.c_int_p = POINTER(c_int)
		self.c_long_p = POINTER(c_long)
		self.c_char_array = ARRAY(c_char,11)
		self.c_char_array_p = POINTER(self.c_char_array)

		# THE FOLLOWING ARE FOR THE PUNCH500 DEL FILES

		# for line parameters
		self.lpars = ([
			'WL',
			'DWL',
			'GFLOG',
			'DGFLOG',
			'CODE',
			'E',
			'XJ',
			'LABEL',
			'EP',
			'XJP',
			'LABELP',
			'GR',
			'DGAMMAR',
			'GS',
			'DGAMMAS',
			'GW',
			'DGAMMAW',
			'WAVENO',
			'REF',
			'NBLO',
			'NBUP',
			'ISO1',
			'X1',
			'ISO2',
			'X2',
			'OTHER1',
			'OTHER2',
			'ISOSHIFT',
			'NELION',
			'RESID',
			])

		# array to convert arrays into actual dtypes
		# keep some arrays as strings because we won't
		# change those.

		self.fmtstr = ({
			'WL':'11.4f',
			'DWL':'7.4f',
			'GFLOG':'7.3f',
			'DGFLOG':'7.3f',
			'CODE':'8.2f',
			'E':'12.3f',
			'XJ':'5.1f',
			'LABEL':'10s',
			'EP':'12.3f',
			'XJP':'5.1f',
			'LABELP':'10s',
			'GR':'6.2f',
			'DGAMMAR':'+6.2f',
			'GS':'6.2f',
			'DGAMMAS':'+6.2f',
			'GW':'6.2f',
			'DGAMMAW':'+6.2f',
			'WAVENO':'11.3f',
			'REF':'4s',
			'NBLO':'2d',
			'NBUP':'2d',
			'ISO1':'3d',
			'X1':'6.3f',
			'ISO2':'3d',
			'X2':'6.3f',
			'OTHER1':'10s',
			'OTHER2':'10s',
			'ISOSHIFT':'6d',
			'NELION':'4d',
			'RESID':'8.4f',
			})
				

		self.fmtnp = {}
		for kk in self.fmtstr.keys():
			if self.fmtstr[kk][-1] == 'f':
				xx = np.float64
			elif self.fmtstr[kk][-1] == 's':
				xx = np.str_
			elif self.fmtstr[kk][-1] == 'd':
				xx = np.int_
			else:
				print('Did not understand format code')
				raise IOError
			self.fmtnp[kk] = xx

		ldelim_i = []
		for kk in self.fmtstr.keys():
			tmpfmt = self.fmtstr[kk][:-1]
			if '.' in tmpfmt:
				num = int(tmpfmt.split('.')[0])
			else:
				num = int(tmpfmt)
			ldelim_i.append(num)
		
		self.ldelim = np.cumsum(ldelim_i)
		self.ldelim_stop = self.ldelim
		self.ldelim_start = [0]
		[self.ldelim_start.append(xx) for xx in self.ldelim]
		self.ldelim_start = self.ldelim_start[:-1]

		self.linedict = {}
		for lp in self.lpars:
			self.linedict[lp] = []



		# THE FOLLOWING ARE FOR THE LINE LIST FROM SYNTOASCANGA.EXE (includes resid at end)

		# for line parameters
		# self.lpars_raw = self.lpars+['RESID']
		
		# array to convert arrays into actual dtypes
		# keep some arrays as strings because we won't
		# change those.
		
		# self.fmtstr_raw = self.fmtstr+['8.4f']
		# self.fmtnp_raw = self.fmtnp+[np.float64]

		# self.ldelim_raw_i = self.ldelim_i+[8]
		
		# self.ldelim_raw = np.cumsum(self.ldelim_raw_i)
		# self.ldelim_raw_stop = self.ldelim_raw
		# self.ldelim_raw_star = [0]
		# [self.ldelim_raw_star.append(xx) for xx in self.ldelim_raw]
		# self.ldelim_raw_star = self.ldelim_raw_star[:-1]

		# self.linedict_raw = {}
		# for lp in self.lpars_raw:
		# 	self.linedict_raw[lp] = []
		
			
	def readspecbin(self,filename,NWL=int(12000),NLINES=int(10000000)):

		s_i  = '{0}'.format(filename)
		s    = bytes(s_i,encoding='ascii')
		SWL   = np.zeros(NWL, dtype="double")
		QMU1 = np.zeros(NWL, dtype="double")
		QMU2 = np.zeros(NWL, dtype="double")

		WLin       = np.zeros(NLINES,dtype='double')
		DWLin      = np.zeros(NLINES,dtype='double')
		GFLOGin    = np.zeros(NLINES,dtype='double')
		DGFLOGin   = np.zeros(NLINES,dtype='double')
		CODEin     = np.zeros(NLINES,dtype='double')
		Ein        = np.zeros(NLINES,dtype='double')
		XJin       = np.zeros(NLINES,dtype='double')
		LABELin    = np.zeros((NLINES,10),dtype='str')
		EPin       = np.zeros(NLINES,dtype='double')
		XJPin      = np.zeros(NLINES,dtype='double')
		LABELPin   = np.zeros((NLINES,10),dtype='str')
		GRin       = np.zeros(NLINES,dtype='double')
		DGAMMARin  = np.zeros(NLINES,dtype='double')
		GSin       = np.zeros(NLINES,dtype='double')
		DGAMMASin  = np.zeros(NLINES,dtype='double')
		GWin       = np.zeros(NLINES,dtype='double')
		DGAMMAWin  = np.zeros(NLINES,dtype='double')
		WAVENOin   = np.zeros(NLINES,dtype='double')
		REFin      = np.zeros((NLINES,4),dtype='str')
		NBLOin     = np.zeros(NLINES,dtype=np.int64)
		NBUPin     = np.zeros(NLINES,dtype=np.int64)
		ISO1in     = np.zeros(NLINES,dtype=np.int64)
		X1in       = np.zeros(NLINES,dtype='double')
		ISO2in     = np.zeros(NLINES,dtype=np.int64)
		X2in       = np.zeros(NLINES,dtype='double')
		OTHER1in   = np.zeros((NLINES,10),dtype='str')
		OTHER2in   = np.zeros((NLINES,10),dtype='str')
		ISOSHIFTin = np.zeros(NLINES,dtype=np.int64)
		NELIONin   = np.zeros(NLINES,dtype=np.int64)
		RESIDin    = np.zeros(NLINES,dtype='double')


		self.fortran.readoutspecbin(
			c_char_p(s), 
			c_int(NWL),
			c_int(NLINES),
			SWL.ctypes.data_as(self.c_double_p),
			QMU1.ctypes.data_as(self.c_double_p),
			QMU2.ctypes.data_as(self.c_double_p),

			WLin.ctypes.data_as(self.c_double_p),      
			DWLin.ctypes.data_as(self.c_double_p),     
			GFLOGin.ctypes.data_as(self.c_double_p),   
			DGFLOGin.ctypes.data_as(self.c_double_p),  
			CODEin.ctypes.data_as(self.c_double_p),    
			Ein.ctypes.data_as(self.c_double_p),       
			XJin.ctypes.data_as(self.c_double_p),      
			LABELin.ctypes.data_as(c_char_p),   
			EPin.ctypes.data_as(self.c_double_p),      
			XJPin.ctypes.data_as(self.c_double_p),     
			LABELPin.ctypes.data_as(c_char_p),  
			GRin.ctypes.data_as(self.c_double_p),      
			DGAMMARin.ctypes.data_as(self.c_double_p), 
			GSin.ctypes.data_as(self.c_double_p),      
			DGAMMASin.ctypes.data_as(self.c_double_p), 
			GWin.ctypes.data_as(self.c_double_p),      
			DGAMMAWin.ctypes.data_as(self.c_double_p), 
			WAVENOin.ctypes.data_as(self.c_double_p),  
			REFin.ctypes.data_as(c_char_p),     
			NBLOin.ctypes.data_as(self.c_long_p),    
			NBUPin.ctypes.data_as(self.c_long_p),    
			ISO1in.ctypes.data_as(self.c_long_p),    
			X1in.ctypes.data_as(self.c_double_p),      
			ISO2in.ctypes.data_as(self.c_long_p),    
			X2in.ctypes.data_as(self.c_double_p),      
			OTHER1in.ctypes.data_as(c_char_p),  
			OTHER2in.ctypes.data_as(c_char_p),  
			ISOSHIFTin.ctypes.data_as(self.c_long_p),
			NELIONin.ctypes.data_as(self.c_long_p),  
			RESIDin.ctypes.data_as(self.c_double_p)  
			)


		x = np.array([''.join(LABELin[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		x = ''.join(x)
		LABELin = list(map(''.join, zip(*[iter(x)]*10)))

		x = np.array([''.join(LABELPin[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		x = ''.join(x)
		LABELPin = list(map(''.join, zip(*[iter(x)]*10)))

		x = np.array([''.join(REFin[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		x = ''.join(x)
		REFin = list(map(''.join, zip(*[iter(x)]*4)))

		x = np.array([''.join(OTHER1in[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		x = ''.join(x)
		OTHER1in = list(map(''.join, zip(*[iter(x)]*10)))

		x = np.array([''.join(OTHER2in[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		x = ''.join(x)
		OTHER2in = list(map(''.join, zip(*[iter(x)]*10)))

		outspec = {'WAVE':SWL,'QMU1':QMU1,'QMU2':QMU2}
		ll = {}
		ll['WL']       = np.array(['{0:11.4f}'.format(x) for x in WLin],dtype='float')
		ll['DWL']      = np.array(['{0:7.4f}'.format(x) for x in DWLin],dtype='float')
		ll['GFLOG']    = np.array(['{0:7.3f}'.format(x) for x in GFLOGin],dtype='float')
		ll['DGFLOG']   = np.array(['{0:7.3f}'.format(x) for x in DGFLOGin],dtype='float')
		# ll['CODE']     = np.array(['{0:8.2f}'.format(x) for x in CODEin],dtype='float')
		ll['CODE']     = np.array(['{0:8.2f}'.format(x) for x in CODEin], dtype='float')
		ll['E']        = np.array(['{0:12.3f}'.format(x) for x in Ein],dtype='float')    
		ll['XJ']       = np.array(['{0:5.1f}'.format(x) for x in XJin],dtype='float')
		# ll['LABEL']    = np.array([''.join(LABELin[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		ll['LABEL']    = np.array([' {0:10s}'.format(x) for x in LABELin],dtype=str)
		ll['EP']       = np.array(['{0:12.3f}'.format(x) for x in EPin],dtype='float')
		ll['XJP']      = np.array(['{0:5.1f}'.format(x) for x in XJPin],dtype='float')
		# ll['LABELP']   = np.array([''.join(LABELPin[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		ll['LABELP']   = np.array([' {0:10s}'.format(x) for x in LABELPin],dtype=str)
		ll['GR']       = np.array(['{0:6.2f}'.format(x) for x in GRin],dtype='float')
		ll['DGAMMAR']  = DGAMMARin
		ll['GS']       = np.array(['{0:6.2f}'.format(x) for x in GSin],dtype='float')
		ll['DGAMMAS']  = DGAMMASin
		ll['GW']       = np.array(['{0:6.2f}'.format(x) for x in GWin],dtype='float')
		ll['DGAMMAW']  = DGAMMAWin
		ll['WAVENO']   = np.array(['{0:11.3f}'.format(x) for x in WAVENOin],dtype='float')
		# ll['REF']      = np.array([''.join(REFin[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		# ll['REF']      = np.array(REFin,dtype='str')
		ll['REF']      = np.array([' {0:4s} '.format(x) for x in REFin],dtype=str)
		ll['NBLO']     = NBLOin
		ll['NBUP']     = NBUPin
		ll['ISO1']     = ISO1in
		ll['X1']       = np.array(['{0:6.3f}'.format(x) for x in X1in],dtype='float')
		ll['ISO2']     = ISO2in
		ll['X2']       = np.array(['{0:6.3f}'.format(x) for x in X2in],dtype='float')  
		# ll['OTHER1']   = np.array([''.join(OTHER1in[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		# ll['OTHER2']   = np.array([''.join(OTHER2in[i,:].tostring('F').decode('ascii')) for i in range(NLINES)])
		ll['OTHER1']   = np.array(OTHER1in,dtype='str')
		ll['OTHER2']   = np.array(OTHER2in,dtype='str')
		ll['ISOSHIFT'] = ISOSHIFTin
		ll['NELION']   = NELIONin
		# ll['OTHER'] = np.array(
		# 	['{0}{1}{2:6.0f}{3:4.0f}'.format(x1,x2,x3,x4) for x1,x2,x3,x4 in zip(
		# 		ll['OTHER1'],ll['OTHER2'],ll['ISOSHIFT'],ll['NELION'])])
		ll['RESID']    = np.array(['{0:8.4f}'.format(x) for x in RESIDin],dtype='float')

		sortind = np.argsort(ll['WL'])

		for kk in ll.keys():
			ll[kk] = ll[kk][sortind]

		for kk in ll.keys():
			for ii in np.random.choice(range(len(ll)),5):
				print(kk,'---','{}'.format(ll[kk][ii]),'---',len('{}'.format(ll[kk][ii])))


		return (outspec,ll)

	def readlp(self,filename):
		"""
		READ(11,140)WL,DWL,GFLOG,DGFLOG,CODE,E,XJ,LABEL,
		1        EP,XJP,LABELP,GR,DGAMMAR,GS,DGAMMAS,GW,DGAMMAW,WAVENO,
		2        REF,NBLO,NBUP,ISO1,X1,ISO2,X2,OTHER1,OTHER2,ISOSHIFT,
		$        NELION
		140     FORMAT(F11.4,F7.4,2F7.3,F8.2,F12.3,F5.1,1X,A8,A2,
		1        F12.3,F5.1,1X,A8,A2,6F6.2,F11.3,
		2        1X,A4,I2,I2,I3,F6.3,I3,F6.3,A8,A2,A8,A2,I6,I4,2X)
		"""

		# make copy of line dictionary
		linedict_i = copy.deepcopy(self.linedict)
		
		# read in ascii line file into an astropy table
		with open(filename,'r') as infile:
			for line in infile:
				for jj,start in enumerate(self.ldelim_star):
					par = line[start:self.ldelim_stop[jj]]
					linedict_i[self.lpars[jj]].append(par)
		# change dictionary to astropy table
		x = Table(linedict_i)
		# order Table correctly because dictionaries are unordered
		x = x[tuple(self.lpars)]
		return x

	def writenp(self,tab,filename):
		# takes an np format table, converts it to lp, and then writes it to file
		tab_i = self.con_nptolp(tab)
		self.writelp(tab_i,filename)
		return None

	def writelp(self,tab,filename):
		# write an ascii line file from an astropy table
		with open(filename,'w') as outfile:
			for ii,tab_i in enumerate(tab):
				if (ii == len(tab)-1):
					outstr = ''.join(tab_i)
				else:
					outstr = ''.join(tab_i)+'\n'					
				outfile.write(outstr)
			outfile.write('\n')
		return None
	
	def con_lptonp(self,strtab):
		outtab = Table(strtab,copy=True,dtype=self.fmtnp)
		return outtab

	def con_nptolp(self,nptab):
		outtab_arr = []
		outpar_arr = []
		for ii,fmt in enumerate(self.lpars):
			if fmt in nptab.dtype.names:
				temparr = ([('{0:'+self.fmtstr[fmt]+'}').format(
							xx) for xx in nptab[fmt]])
				outtab_arr.append(temparr)
				outpar_arr.append(fmt)
		outtab = Table(outtab_arr,names=outpar_arr)
		return outtab

	def readlp_raw(self,filename):
		"""
		read in raw line list (i.e., no delta values). Takes in 
		ascii file from syntoascanga.exe
		"""
		# make a copy of line dictionary
		self.linedict_raw_i = copy.deepcopy(self.linedict_raw)

		# read in ascii raw line list into astropy table
		with open(filename,'r') as infile:
			for line in infile:
				for jj,start in enumerate(self.ldelim_raw_star):
					par = line[start:self.ldelim_raw_stop[jj]]
					self.linedict_raw_i[self.lpars_raw[jj]].append(par)

		# change dictionary to astropy table
		x = Table(self.linedict_raw_i)
		# order Table correctly because dictionaries are unordered
		x = x[tuple(self.lpars_raw)]

		return x

	def con_lptonp_raw(self,strtab):
		outtab = Table(strtab,copy=True,dtype=self.fmtnp_raw)
		return outtab
		
