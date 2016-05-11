import numpy as np

# FUNCTION THAT TAKES A FM.LL + CONDITIONAL STATEMENT AND RETURNS A P ARRAY, P_SIGMA ARRAY, AND T ARRAY

def linesel(LL,condst,minWL,maxWL):

	# GENERATE ZEROED T ARRAY
	t = np.array([[-1,-1,-1,-1,-1] for _ in LL])

	# initialize p and p_sig array
	p = []
	psig = []
	pflag = []

	# parameter sigmas to use
	# init_dellam_sig = 0.01
	# init_delgf_sig  = 2.5
	# init_delga_sig = 2.5

	# FOR TESTING PURPOSES
	init_dellam_sig = 0.01
	init_delgf_sig  = 1.5
	init_delga_sig = 0.5

	# find the lines to be fit
	if condst == None:
		# using all of the lines within this line segment
		cond_tot = np.less(LL['WL'],maxWL) * np.greater(LL['WL'],minWL)
		cond_tot = np.array(cond_tot)
		lineselectind = np.argwhere(cond_tot).ravel()

	else:
		# using conditional string
		cond_tot = [True for _ in LL['WL']]
		for condict in condst:
			cond_i = condict['OP'](LL[condict['LP']],condict['LV'])
		cond_tot = cond_i * cond_tot
		cond_tot = cond_tot * np.less(LL['WL'],maxWL) * np.greater(LL['WL'],minWL)
		cond_tot = np.array(cond_tot)

		lineselectind = np.argwhere(cond_tot).ravel()

	kk = 0

	# generate band information for mol lines
	molbandarr = np.empty(len(LL),dtype=object)

	for ii,code in enumerate(LL['CODE']):
		if float(code) > 100.0:
			# pick off H2O lines
			if float(code) = 10108.0:
				molbandarr[ii] = 'H2O'
			else:
				# find first character in label 
				for aa in LL['LABEL'][ii]:
					if aa.isalpha():
						AA = aa
						break
				# find first character in labelp 
				for bb in LL['LABELP'][ii]:
					if bb.isalpha():
						BB = bb
						break
				# generate a band symbol
				try:
					molbandarr[ii] = AA+BB
				except UnboundLocalError:
					print LL[ii]
					raise


	for ii,lind in enumerate(lineselectind):
		# CHECK TO MAKE SURE LINE HAS NOT ALREADY BEEN SET BECAUSE IT IS A MOLE/HF/ISO OF PREVIOUS LINE
		if all(t[lind] == -1):
			# append a WL shift
			p.append(LL['DWL'][lind])
			psig.append(init_dellam_sig)
			pflag.append('WL')
			t[lind][0] = kk
			kk = kk+1

			# append a log(gf) shift
			p.append(LL['DGFLOG'][lind])
			psig.append(init_delgf_sig)
			pflag.append('GF')
			t[lind][1] = kk
			kk = kk+1

			# append van der Waals if it isn't a molecule
			if float(LL['CODE'][lind]) < 100.0:
				# van der Waals
				p.append(LL['DGAMMAW'][lind])
				psig.append(init_delga_sig)
				pflag.append('GW')
				t[lind][2] = kk
				kk = kk + 1

			# if fm.gammaswitch['switch'][lind] == 'W':
			# 	# van der Waals
			# 	p.append(LL['DGAMMAW'][lind])
			# 	psig.append(init_delga_sig)
			# 	pflag.append('GW')
			# 	t[lind][2] = kk
			# 	kk = kk + 1

			# elif fm.gammaswitch['switch'][lind] == 'R':
			# 	# radiative
			# 	p.append(LL['DGAMMAR'][lind])
			# 	psig.append(init_delga_sig)
			# 	pflag.append('GR')
			# 	t[lind][3] = kk
			# 	kk = kk + 1

			# elif fm.gammaswitch['switch'][lind] == 'S':
			# 	# starks
			# 	p.append(LL['DGAMMAS'][lind])
			# 	psig.append(init_delga_sig)
			# 	pflag.append('GS')
			# 	t[lind][4] = kk
			# 	kk = kk + 1

			# else:
			# 	pass

			# IDENTIFY ANY POSSIBLE COUPLED LINES AND SET THEIR T ARRAY TO BE THE SAME AS THE CURRENT LINE

			# FIRST SEE IF LIND IS A MOLECULE USING GAMMASWITCH AND SET ALL MOLECULES OF SAME TYPE IN WAVERANG TO SAME LOG(gf)
			if float(LL['CODE'][lind]) > 100.0:
				molind = np.where( 
					(LL['CODE'] == LL['CODE'][lind]) & 
					(LL['ISO1'] == LL['ISO1'][lind]) &
					(LL['ISO2'] == LL['ISO2'][lind]) &
					(molbandarr == molbandarr[lind]) &
					np.less(LL['WL'],maxWL) & np.greater(LL['WL'],minWL)
					)
				t[molind] = t[lind]
				# CONTINUE SINCE NO ISOTOPIC OR HFS OF MOLECULES
				continue
			# SECOND IDENTIFY ANY POSSIBLE HF OR ISO LINES BASED ON HAVING THE SAME E, EP, LABEL, and LABELP
			HFISOind = np.where( 
				(LL['E'] == LL['E'][lind]) & 
				(LL['EP'] == LL['EP'][lind]) & 
				(LL['LABEL'] == LL['LABEL'][lind]) & 
				(LL['LABELP'] == LL['LABELP'][lind]) &
				np.less(LL['WL'],maxWL) & np.greater(LL['WL'],minWL) 
				)
			# set their pars to the current lines
			if len(HFISOind[0]) > 0:
				t[HFISOind] = t[lind]

	return (np.array(p),np.array(psig),np.array(pflag),t)