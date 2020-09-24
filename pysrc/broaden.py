from ctypes import cdll, CDLL, POINTER, c_int, c_double, c_char_p, c_char, c_void_p
import numpy as np
import os

import socket
hostname = socket.gethostname()
if hostname[:4] == 'holy':
     RUNLOC = 'ODY'
else:
     RUNLOC = 'LOCAL'

conroypath = os.environ.get('CSCRATCH')
holypath   = os.environ.get('HOLYSCRATCH')
homepath   = os.environ.get('HOME')

class broaden(object):
	def __init__(self):
		"""
		Class for python hooks for Bob Kurucz broadening function.
		"""

		# read in fortran library with ctype hooks
		if RUNLOC == 'ODY':
			self.fortran = cdll.LoadLibrary('{0}/pac/broadenx-py/libbroaden.so'.format(conroypath))		
		else:
			self.fortran = cdll.LoadLibrary('/Users/pcargile/Astro/bin/FORTRAN/SYNTHE/lib/libbroaden.so')

		# define this useful thing for later
		self.c_double_p = POINTER(c_double)

	def broaden(self,wli,fluxi,broaddict):
		"""
		Function to take in wavelength and flux arrays 
		and output broaden arrays using the fortran 
		broadening code.

		wli = input wavelength array
		fluxi = input flux array
		resol = resolution of input spectrum
		broaddict = dictionary of info for broadening
				possible dictionary keys:
				type: (required) type of broadening (MACRO, GAUSSIAN, SINX/X, RECT)
				units: (required) units for broadening:
				 	PM
				 	KM
				 	RESOLUTION
				 	CM-1
				and either
				val: broadening over full spectrum
				or using a linear change in resolution or resolving power over spectrum
				val1: left edge broadening value
				val2: right edge broadening value

		OUTPUT:
			dictionary with "WAVE" and "FLUX" broadened spectrum

		"""

		# make sure input arrays are np arrays and dtype=double
		wli = np.array(wli,dtype="double")
		fluxi = np.array(fluxi,dtype="double")

		# define other values used in fortan call
		N = len(wli) # number of rows in wave/flux arrays
		A = bytes('{0:10s}'.format( broaddict['type']),encoding='ascii') # type of broadening
		B = bytes('{0:10s}'.format(broaddict['units']),encoding='ascii') # units on broadening type
		if 'val' in broaddict.keys():
			X1 = float(broaddict['val']) # single broadening value
			X2 = float(0.0) # set second to zero to turn off
		if 'val1' in broaddict.keys():
			X1 = float(broaddict['val1']) # first broadening value
		if 'val2' in broaddict.keys():
			X2 = float(broaddict['val2']) # second broadening value
		else:
			X2 = float(0.0) # in case there is only one broadening value

		# define output array
		fluxo = np.zeros(N,dtype="double")

		# run the fortran broadening code
		self.fortran.broaden(
			c_char_p(A), 
			c_double(X1), 
			c_char_p(B),
			c_double(X2),
			c_int(N), 
			# c_double(resol), 
			wli.ctypes.data_as(self.c_double_p), 
			fluxi.ctypes.data_as(self.c_double_p), 
			fluxo.ctypes.data_as(self.c_double_p)
			)
		# return dictionary with broadened flux values
		return {"WAVE":wli,"FLUX":fluxo}
