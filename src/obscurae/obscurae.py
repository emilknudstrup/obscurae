#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# =============================================================================
# Packages
# =============================================================================

import numpy as np
import scipy.signal as ss

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# =============================================================================
# Grid
# =============================================================================

class Disk(object):
	'''Disk grid.

	Instantiate the stellar disk.

	:param npix: Number of pixels in the stellar disk. Default is 201.
	:type npix: int, optional
	
	:param thick: Thickness of rings. Default is 20 pixels).
	:type thick: int, optional

	:param span: Span of velocity grid. Default is 20 km/s.
	:type span: int, optional

	:param res: Resolution of velocity grid. Default is 0.25 km/s.
	:type res: float, optional

	'''


	def __init__(self,
		npix = 201,
		thick = 20,
		span = 20,
		res = 0.25
		):
		self.npix = npix
		self.thick = thick
		self.vel_res = np.arange(-span,span,res)

	@staticmethod
	def gridCoords(rad,xoff=0,yoff=0):
		'''Coordinates of the grid.
		
		Calculate the coordinates of the grid, from -rad to +rad.

		:param rad: Radius of the grid in number of pixels.
		:type rad: int

		:param xoff: Offset in x-direction. Default is 0.
		:type xoff: float, optional

		:param yoff: Offset in y-direction. Default is 0.
		:type yoff: float, optional

		:return: grid coordinates, indices of coordinates, distance from limb
		:rtype: (array, array, array)

		'''

		xx, yy = np.arange(-rad+xoff,rad+1+xoff), np.arange(-rad+yoff,rad+1+yoff)
		coord = np.array(np.meshgrid(xx,yy)).T.reshape(-1,2)

		rad_dist = np.sqrt(np.add.reduce(coord**2,axis=1))

		dd = np.arange(0,2*rad+1,1)
		cidxs = [tuple(cc) for cc in np.array(np.meshgrid(dd,dd)).T.reshape(-1,2)] # Indices for coordinates

		return coord, cidxs, rad_dist

	@staticmethod
	def grid(rad,xoff=0,yoff=0):
		'''Initial grid.

		:param rad: Radius in number of pixels.
		:type rad: int 
		
		:param xoff: Offset in x-direction. Default 0.
		:type xoff: float, optional 
		
		:param yoff: Offset in y-direction. Default 0.
		:type yoff: float, optional 

		:return: initial grid, velocity grid, normalized radial coordinate.
		:rtype: (array, array, array)

		'''
		
		coord, cidxs, rad_dist = Disk.gridCoords(rad,xoff=xoff,yoff=yoff)
		## Empty, quadratic array with dimensions diameter*diameter
		start_grid = np.zeros((2*rad+1,2*rad+1))
		vel = start_grid.copy()
		mu = start_grid.copy()
		for ii in range(len(cidxs)):
			if rad_dist[ii] <= rad:
				start_grid[cidxs[ii]] = 1
				mu[cidxs[ii]] = np.sqrt(1-(rad_dist[ii]/float(rad))**2)
			vel[cidxs[ii]] = coord[ii][0]/float(rad)
		
		return start_grid, vel, mu

	@staticmethod
	def gridRing(rad,thick,xoff=0,yoff=0):
		'''Initial grid in rings of :math:`\mu`.

		Initial grid of star in rings of aprrox same :math:`\mu = \cos(\\theta)`, 
		where :math:`\\theta` is the angle between a line through the center of the star and the limb.

		Useful for macroturbulence calculations.
		
		:param Rs: Radius in number of pixels.
		:type Rs: int 
		
		:param thick: Thickness of rings.
		:type thick: int
		
		:param xoff: Potential offset in x-direction. Default 0.
		:type xoff: float, optional 
		
		:param yoff: Potential offset in y-direction. Default 0.
		:type yoff: float, optional 

		:return: pixels within stellar disk, velocity grid, radial :math:`\mu` values, approx :math:`\mu` in each ring
		:rtype: (array, array, array, array)

		'''

		assert rad/thick > 2.0, print('The radius must be at least twice the size of the rings.')
		coord, cidxs, rad_dist = Disk.gridCoords(rad,xoff=xoff,yoff=yoff)

		## Empty, quadratic array with dimensions diameter*diameter
		start_grid = np.zeros((2*rad+1,2*rad+1))
		vel = start_grid.copy()	

		## Divide star into rings
		rings = np.arange(0,rad+1,thick)
		nr = len(rings) # number of rings
		rings_in = rings - thick
		## Make sure to get all of the star
		if rings[-1] < rad:
			rings_in = np.append(rings_in,rings[-1])
			rings = np.append(rings,rad)

		mu_grid = np.asarray([start_grid.copy() for i in range(nr)])
		ring_grid = mu_grid.copy()#np.asarray([start_grid.copy() for i in range(nr)]) #startgrid for ring
		for ii in range(len(cidxs)):
			for jj in range(nr):
				if rings_in[jj] < rad_dist[ii] <= rings[jj]:
					ring_grid[jj][cidxs[ii]] = 1
					mu_grid[jj][cidxs[ii]] = np.sqrt(1-(rad_dist[ii]/float(rad))**2) ## Used for limb darkening and macro. Same as cos(theta)
			vel[cidxs[ii]] = coord[ii][0]/float(rad) ## Velocity field. [0]-axis (rows) are cst vel, while [1]-axis (columns) are the equator
		
		mu_mean = np.zeros(shape=(1,nr))[0]
		## Calculate the approx mu in each ring to be used when calculating the macroturbulence
		for kk in range(nr):
			mu_mean[kk] = np.mean(mu_grid[kk][np.where(mu_grid[kk])])


		return ring_grid, vel, mu_grid, mu_mean




	def Grid(self):
		'''Instantiate the grid.

		

		'''
		

		## Initial grid
		self.star_grid, _, self.mu = self.grid(self.npix)	
		## The grid is made into rings for faster calculation of the macroturbulence (approx. constant in each ring)
		self.ring_grid, self.vel, self.mu_grid, self.mu_mean = self.gridRing(self.npix,self.thick)


	def showGrid(self,ax=None,raw=False,rings=False,ring_cmap='copper'):
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		if raw:
			ax.imshow(self.star_grid.T, cmap='Greys_r', interpolation='nearest')
		else:
			ax.imshow(self.mu.T, cmap='Greys_r', interpolation='nearest')	
		if rings:
			nr = self.ring_grid.shape[0]
			cmap = plt.get_cmap(ring_cmap,nr)
			for ii in range(nr):
				ax.imshow(self.ring_grid[ii],cmap=ListedColormap(['none', cmap(ii)]),alpha=0.35)
	
		ax.imshow(self.star_grid.T, cmap=ListedColormap(['k', 'none']), interpolation='nearest')

		ax.set_ylim(ax.get_ylim()[::-1])
		ax.set_xlabel('x (pixels)')
		ax.set_ylabel('y (pixels)')

# =============================================================================
# Stellar grid
# =============================================================================


class Star(Disk):
	'''Stellar grid.

	Instantiate the stellar grid.
	Inherits from :class:`Disk`.

	:param vsini: Stellar rotation velocity. Default is 10 km/s.
	:type vsini: float, optional

	:param zeta: Macroturbulence parameter. Default is 2 km/s.
	:type zeta: float, optional

	:param xi: Microturbulence parameter. Default is 1 km/s.
	:type xi: float, optional

	:param cs: Limb darkening coefficients. Default is [0.6,0.4].
	:type cs: list, float, optional

	:param law: Limb darkening law. Default is quadratic.
	:type law: str, optional

	:param npix: Number of pixels in the stellar disk. Default is 201.
	:type npix: int, optional

	:param thick: Thickness of rings. Default is 20 pixels).
	:type thick: int, optional

	:param span: Span of velocity grid. Default is 20 km/s.

	:param res: Resolution of velocity grid. Default is 0.25 km/s.
	:type res: float, optional

	'''


	def __init__(self,
		vsini=10.0,zeta=2.0,xi=1.0,cs=[0.6,0.4],law='quad',#star
		npix = 100,thick = 20,span = 20,res = 0.25#grid
		):

		self.vsini = vsini
		self.zeta = zeta
		self.xi = xi
		self.cs = cs
		self.law = law

		Disk.__init__(self,npix,thick,span,res)


	def limbDarkening(self,cs=[0.3,0.2],law='quad'):
		'''Limb darkening.

		Calculates the limb darkening of the star based on :math:`\mu` and the limb darkening law.

		:param cs: Limb darkening coefficients. Default is [0.3,0.2].
		:type cs: list, float

		:param law: Limb darkening law. Default is quadratic.
		:type law: str, optional

		.. note::
			The 'nonlinear' law needs four coefficients, while the 'quad' only need two. 'uni' naturally does not need any.

		'''
		

		if law == 'quad':
			c1, c2 = cs
			LD = 1 - c1*(1 - self.mu) - c2*np.power(1 - self.mu,2)
		elif law == 'nonlinear':
			c1, c2, c3, c4 = cs
			LD = 1. - c1*(1. - np.sqrt(self.mu)) - c2*(1. - self.mu) - c3*(1. - np.power(self.mu,3/2)) - c4*(1. - np.power(self.mu,2))
		elif law == 'uni':
			LD = 1

		self.LD_grid = self.ring_grid*LD


	def microturbulence(self,xi):
		'''Microturbulence.

		Convolves the rotation profile with a gaussian to take microturbulence for a given :math:`\\xi` into account.
		
		:param xi: Microturbulence parameter.
		:type xi: float

		'''

		## Gaussian function for microturbuelence
		## make Gaussian with new velocity vector as x-axis
		gau = np.exp(-1*np.power(self.vel_res/xi,2))/(xi*np.sqrt(np.pi))#gauss(x,xi) 
		gau /= np.add.reduce(gau)

		length = len(self.active_ring[0])+len(gau)-1
		gau_arr = np.empty([len(self.active_ring),len(gau)])
		gau_arr[:] = gau

		line = np.add.reduce(self.active_ring,axis=2)
		micro_profile = np.empty([len(self.active_ring),length])
		micro_profile[:] = ss.fftconvolve(line[:],gau_arr[:],'full',axes=1)
		self.micro = micro_profile


	def macroturbulence(self,zeta):
		'''Macroturbulence.

		Calculates the macroturbulence of the star in rings 
		at a given :math:`\mu` 
		for a given :math:`\zeta` using the radial-tangential profile.

		:param zeta: Macroturbulence parameter.
		:type zeta: float

		'''

		## 1D-velocity (same for each ring)
		vel_1d = self.rot_profile[:,0]


		## Calculate macroturbulence of rings
		A = 0.5 #area of curve covered by radial and tangential
		macro_profile = np.zeros(shape=(len(self.mu_mean),len(vel_1d)))	

		for ii, mu in enumerate(self.mu_mean):
			rad = np.exp(-1*np.power(vel_1d/(zeta*mu),2))/mu
			rad[np.isnan(rad)] = 0.
			y = np.sin(np.arccos(mu))
			tan = np.exp(-1*np.power(vel_1d/(zeta*y),2))/y
			tan[np.isnan(tan)] = 0.
			macro_profile[ii] = A*(rad + tan)/(np.sqrt(np.pi)*zeta)

		self.macro = macro_profile

	def convolve(self):
		'''Convolve the stellar grid.
		
		Convolces the stellar grid with the macroturbulence and microturbulence following the approaches in :cite:t:`Hirano2011` and :cite:t:`Gray2005`.
		
		'''

		## 1D-velocity (same for each ring)
		vel_1d = self.rot_profile[:,0]
		sep = (vel_1d[-1]-vel_1d[0])/(len(vel_1d)-1)

		## Convolve each ring of LD+microturbulence and macroturbulence
		lc_len = self.macro[0].shape[0] + self.micro[0].shape[0] - 1 
		xn = np.linspace(-lc_len*sep/2.,lc_len*sep/2.,num=lc_len,endpoint=False)

		line_conv = np.zeros(shape=(len(self.active_ring),lc_len))
		line_conv[:,:] = ss.fftconvolve(self.micro,self.macro,'full',axes=1)

		self.vel_1d = xn
		self.line_conv = line_conv


	def Line(self):
		'''Stellar line profile.

		The stellar line profile is calculated by convolving the
		rotating, limb-darkened star grid with the macroturbulence
		and microturbulence.


		'''

		## Rotation profile
		self.rot_profile = self.vel*self.vsini

		## Limb darkened disk
		self.limbDarkening(cs=self.cs,law=self.law)

		## Make the limb-darkened (LD) stellar grid into rings
		LD_ring = self.ring_grid*self.LD_grid
		self.active_ring = LD_ring.copy()

		## Microturbulence
		self.microturbulence(self.xi)

		## Macroturbulence
		self.macroturbulence(self.zeta)

		## Stellar line profile
		self.convolve()

		line_oot = np.sum(self.line_conv,axis=0)
		self.area = np.trapz(line_oot,self.vel_1d)
		self.line_profile = line_oot/self.area


	## plot the star
	def starLight(self,ax=None,raw=False,rotation=True):
		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)

		if rotation:
			rot = np.multiply(self.star_grid,self.rot_profile)#,out=self.star_grid)
			ax.imshow(rot.T, cmap='coolwarm', interpolation='nearest')	
			ax.imshow(np.add.reduce(self.LD_grid,axis=0).T, cmap='Greys_r', interpolation='nearest',alpha=0.5)
		else:
			ax.imshow(np.add.reduce(self.LD_grid,axis=0).T, cmap='afmhot', interpolation='nearest',alpha=0.75)


		ax.imshow(self.star_grid.T, cmap=ListedColormap(['k', 'none']), interpolation='nearest')

		ax.set_ylim(ax.get_ylim()[::-1])
		ax.set_xlabel('x (pixels)')
		ax.set_ylabel('y (pixels)')

	def starLine(self,ax=None):

		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)

		ax.plot(self.vel_1d,self.line_profile,color='k',lw=3.0)
		ax.plot(self.vel_1d,self.line_profile,color='C0',lw=2.0)
		ax.set_xlabel('Velocity (km s$^{-1}$)')
		ax.set_ylabel('CCF')


# =============================================================================
# Planet grid
# =============================================================================


class Planet(Star):
	'''Planet grid.

	Instantiate the planet grid.
	Inherits from :class:`Star`.

	:param rp: Planet radius in stellar radii. Default is 0.1.
	:type rp: float, optional

	:param vsini: Stellar rotation velocity. Default is 10 km/s.
	:type vsini: float, optional

	:param zeta: Macroturbulence parameter. Default is 2 km/s.
	:type zeta: float, optional

	:param xi: Microturbulence parameter. Default is 1 km/s.
	:type xi: float, optional

	:param cs: Limb darkening coefficients. Default is [0.6,0.4].
	:type cs: list, float, optional

	:param law: Limb darkening law. Default is quadratic.
	:type law: str, optional

	:param npix: Number of pixels in the stellar disk. Default is 201.
	:type npix: int, optional

	:param thick: Thickness of rings. Default is 20 pixels).
	:type thick: int, optional

	:param span: Span of velocity grid. Default is 20 km/s.
	:type span: int, optional

	:param res: Resolution of velocity grid. Default is 0.25 km/s.
	:type res: float, optional

	'''

	def __init__(self,
			rp=0.1,#planet
			vsini=10.0, zeta=2.0, xi=1.0, cs=[0.6,0.4], law='quad',#star
			npix=201, thick=20, span=20, res=0.25#grid
			):
		
		self.rp = rp
		self.is_transiting = False

		Star.__init__(self,
				vsini=vsini,zeta=zeta,xi=xi,cs=cs,law=law,
				npix=npix,thick=thick,span=span,res=res)



	def Occult(self,xx,yy):
		'''Transiting planet on stellar disk.

		The limb-darkened stellar disk is occulted by a planet.
		We grab those pixels occulted by the planet and make a new (sub)grid.

		:param xx: x-position of planet on stellar disk.
		:type xx: array, float

		:param yy: y-position of planet on stellar disk.
		:type yy: array, float

		'''
		if type(xx) == float:
			xx, yy = [xx], [yy]
		xx, yy = np.asarray(xx), np.asarray(yy)
		## grid of planet
		rad = int(round(self.rp*self.npix))
		pl_grid, pl_vel, mu = self.grid(rad)
		rescape_pl = np.reshape(pl_grid,np.size(pl_grid))

		xnorm, ynorm = xx*self.npix, yy*self.npix
		x_off, y_off = np.rint(xnorm), np.rint(ynorm)
		x_pl, y_pl = abs(x_off)-rad, abs(y_off)-rad

		nn = len(xx)
		self.transits = [self.star_grid.copy() for i in range(len(xx))]
		self.path = self.star_grid.copy()
		
		self.planet_rings = [np.zeros_like(self.LD_grid) for i in range(nn)]
		
		for ii in range(nn):
			if (x_pl[ii] < self.npix) & (y_pl[ii] < self.npix):
				
				x_pos, y_pos = int(x_off[ii]+self.npix), int(y_off[ii]+self.npix)

				pl_coord, pl_coord_arr, pl_coord_idx = Disk.gridCoords(rad,x_pos,y_pos)
				
				## Makes sure that when a planet is about to disappear 
				## it does not appear on the other side.
				pl_zip = [(pl_coord[ii,0],pl_coord[ii,1]) for ii in range(pl_coord.shape[0])]
				coord_pl = [i for (i,v) in zip(pl_zip,rescape_pl) if v==1. and i[0] >= 0 and i[1] >= 0 and i[0] < self.LD_grid[0].shape[0] and i[1] < self.LD_grid[0].shape[0]]
				
				self.planet_rings[ii][:,np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]] = self.LD_grid[:,np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]]

				## Only used for plotting
				self.path[np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]] = 0
				self.transits[ii][np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]] = 0.#self.LD_grid[np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]]
				self.is_transiting = True
				

	def transit(self,ax=None,rotation=True,trace=()):

		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		if rotation:
			rot = np.multiply(self.star_grid,self.rot_profile)#,out=self.star_grid)
			ax.imshow(rot.T, cmap='coolwarm', interpolation='nearest')	
			ax.imshow(np.add.reduce(self.LD_grid,axis=0).T, cmap='Greys_r', interpolation='nearest',alpha=0.5)
		else:
			ax.imshow(np.add.reduce(self.LD_grid,axis=0).T, cmap='afmhot', interpolation='nearest',alpha=0.75)

		if self.is_transiting:
			ax.imshow(self.path.T, cmap=ListedColormap(['k', 'none']), interpolation='nearest')
		
		ax.set_ylim(ax.get_ylim()[::-1])
		ax.set_xlim(ax.get_xlim())
		if len(trace):
			xx, yy = trace
			ax.plot(xx*self.npix+self.npix,yy*self.npix+self.npix,'k-',lw=2.0)
			ax.plot(xx*self.npix+self.npix,yy*self.npix+self.npix,'w-',lw=1.0)
		ax.set_xlabel('x (pixels)')
		ax.set_ylabel('y (pixels)')		

		ax.imshow(self.star_grid.T, cmap=ListedColormap(['k', 'none']), interpolation='nearest')


	def distortLine(self):

		## Rotation profile
		self.rot_profile = self.vel*self.vsini

		## Limb darkened disk
		self.limbDarkening(cs=self.cs,law=self.law)

		## Make the limb-darkened (LD) stellar grid into rings
		LD_ring = self.ring_grid*self.LD_grid
		self.active_ring = LD_ring.copy()

		## maybe not needed
		self.active_mu = self.mu_mean.copy()

		## Macroturbulence
		self.macroturbulence(self.zeta)

		## Microturbulence
		self.microturbulence(self.xi)

		## Stellar line profile
		self.convolve()

		line_oot = np.sum(self.line_conv,axis=0)
		self.area = np.trapz(line_oot,self.vel_1d)
		self.line_profile = line_oot/self.area

		## Distorted line profiles
		self.distorted_lines = []
		for ii in range(len(self.planet_rings)):

			## Make the limb-darkened (LD) stellar grid into rings
			LD_ring = self.planet_rings[ii]*self.ring_grid
			self.active_ring = LD_ring.copy()

			## Microturbulence
			## only micro-turbulence is different
			## WHEN implemented this way
			self.microturbulence(self.xi)

			## Stellar line profile
			self.convolve()

			line = np.add.reduce(self.line_conv,axis=0)
			dline = line_oot - line
			self.distorted_lines.append(dline/self.area)

	def showLines(self,ax=None):

		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)

		for ii in range(len(self.distorted_lines)):
			ax.plot(self.vel_1d,self.line_profile-self.distorted_lines[ii],color='k',lw=3.0)
			ax.plot(self.vel_1d,self.line_profile-self.distorted_lines[ii],lw=2.0)
		ax.set_xlabel('Velocity (km s$^{-1}$)')
		ax.set_ylabel('CCF')




# =============================================================================
# Spot
# =============================================================================

class Spot(Star):
	'''Spot grid.

	Instantiate the spot grid.
	Inherits from :class:`Star`.

	:param Tspot: Spot temperature in K. Default is 2500 K.
	:type Tspot: float, optional

	:param Teff: Stellar effective temperature in K. Default is 5777 K.
	:type Teff: float, optional

	:param Rspot: Spot radius in stellar radii. Default is 0.1.
	:type Rspot: float, optional

	:param vsini: Stellar rotation velocity. Default is 10 km/s.
	:type vsini: float, optional

	:param zeta: Macroturbulence parameter. Default is 2 km/s.
	:type zeta: float, optional

	:param xi: Microturbulence parameter. Default is 1 km/s.
	:type xi: float, optional

	:param cs: Limb darkening coefficients. Default is [0.6,0.4].
	:type cs: list, float, optional

	:param law: Limb darkening law. Default is quadratic.
	:type law: str, optional

	:param npix: Number of pixels in the stellar disk. Default is 201.
	:type npix: int, optional

	:param thick: Thickness of rings. Default is 20 pixels).
	:type thick: int, optional

	:param span: Span of velocity grid. Default is 20 km/s.
	:type span: int, optional

	:param res: Resolution of velocity grid. Default is 0.25 km/s.
	:type res: float, optional

	'''

	def __init__(self,
		Tspot=2500, Teff=5777, Rspot=0.1, #spot
		vsini=10.0,zeta=2.0,xi=1.0,cs=[0.6,0.4],law='quad',#star
		npix = 100,thick = 20,span = 20,res = 0.25#grid
		):

		self.Tspot = Tspot
		self.Teff = Teff
		self.Rspot = Rspot

		Star.__init__(self,
				vsini=vsini,zeta=zeta,xi=xi,cs=cs,law=law,
				npix=npix,thick=thick,span=span,res=res)


	def Cross(self,xx,yy):
		'''Spot crossing stellar disk.

		Similar to :func:`Planet.Occult` but for spots with fractional intensity.
		
		:param xx: x-position of spot on stellar disk.
		:type xx: array, float

		:param yy: y-position of spot on stellar disk.
		:type yy: array, float

		'''
		## Check if xx and yy are floats,
		## convert to arrays
		if type(xx) == float:
			xx, yy = [xx], [yy]
		xx, yy = np.asarray(xx), np.asarray(yy)
		## grid of planet
		rad = int(round(self.Rspot*self.npix))
		sp_grid, sp_vel, mu = self.grid(rad)
		rescape_sp = np.reshape(sp_grid,np.size(sp_grid))

		## Fraction of flux blocked by spot
		frac = np.power(self.Tspot/self.Teff,4)

		xnorm, ynorm = xx*self.npix, yy*self.npix
		x_off, y_off = np.rint(xnorm), np.rint(ynorm)
		x_pl, y_pl = abs(x_off)-rad, abs(y_off)-rad

		nn = len(xx)

		self.spots = [self.star_grid.copy() for i in range(nn)]
		self.path = self.star_grid.copy()
		self.spot_rings = [np.zeros_like(self.LD_grid) for i in range(nn)]
		
		for ii in range(nn):
			if (x_pl[ii] < self.npix) & (y_pl[ii] < self.npix):
				
				x_pos, y_pos = int(x_off[ii]+self.npix), int(y_off[ii]+self.npix)
				sp_coord, _, _ = Disk.gridCoords(rad,x_pos,y_pos)
	
				## Makes sure that when a spot is about to disappear 
				## it does not appear on the other side.
				pl_zip = [(sp_coord[ii,0],sp_coord[ii,1]) for ii in range(sp_coord.shape[0])]
				coord_pl = [i for (i,v) in zip(pl_zip,rescape_sp) if v==1. and i[0] >= 0 and i[1] >= 0 and i[0] < self.LD_grid[0].shape[0] and i[1] < self.LD_grid[0].shape[0]]
				
				self.spot_rings[ii][:,np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]] = self.LD_grid[:,np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]]*(1-frac)

				self.spots[ii][np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]] = frac
				self.path[np.asarray(coord_pl)[:,0],np.asarray(coord_pl)[:,1]] = frac		


	def crossing(self,ax=None,rotation=True,trace=()):

		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)
		if rotation:
			rot = np.multiply(self.star_grid,self.rot_profile)#,out=self.star_grid)
			ax.imshow(rot.T, cmap='coolwarm', interpolation='nearest')	
			ax.imshow(np.add.reduce(self.LD_grid,axis=0).T, cmap='Greys_r', interpolation='nearest',alpha=0.5)
		else:
			ax.imshow(np.add.reduce(self.LD_grid,axis=0).T, cmap='afmhot', interpolation='nearest',alpha=0.75)

		#ax.imshow(self.spots.T, cmap=ListedColormap(['k','none']), interpolation='nearest',alpha=np.power(self.Tspot/self.Teff,4))

		#ax.imshow(self.path.T, cmap=ListedColormap(['k','dimgray','darkgray','lightgray','none']), interpolation='nearest',alpha=0.9)
		#frac = np.power(self.Tspot/self.Teff,4)
		frac = self.Tspot/self.Teff
		if frac < 1.0:
			ax.imshow(self.path.T, cmap=ListedColormap(['C3','none']), interpolation='nearest',alpha=1-frac)
		elif frac < 2.0:
			ax.imshow(self.path.T, cmap=ListedColormap(['none','C0']), interpolation='nearest',alpha=frac-1.0)
		else:
			ax.imshow(self.path.T, cmap=ListedColormap(['none','w']), interpolation='nearest',alpha=1.0)

		ax.set_ylim(ax.get_ylim()[::-1])
		ax.set_xlim(ax.get_xlim())
		if len(trace):
			xx, yy = trace
			ax.plot(xx*self.npix+self.npix,yy*self.npix+self.npix,'k-',lw=2.0)
			ax.plot(xx*self.npix+self.npix,yy*self.npix+self.npix,'w-',lw=1.0)

		ax.imshow(self.star_grid.T, cmap=ListedColormap(['k', 'none']), interpolation='nearest')

		ax.set_xlabel('x (pixels)')
		ax.set_ylabel('y (pixels)')

	def distortLine(self):

		## Rotation profile
		self.rot_profile = self.vel*self.vsini

		## Limb darkened disk
		self.limbDarkening(cs=self.cs,law=self.law)

		## Make the limb-darkened (LD) stellar grid into rings
		LD_ring = self.ring_grid*self.LD_grid
		self.active_ring = LD_ring.copy()

		## maybe not needed
		self.active_mu = self.mu_mean.copy()

		## Macroturbulence
		self.macroturbulence(self.zeta)

		## Microturbulence
		self.microturbulence(self.xi)

		## Stellar line profile
		self.convolve()

		line_oot = np.sum(self.line_conv,axis=0)
		self.area = np.trapz(line_oot,self.vel_1d)
		self.line_profile = line_oot/self.area

		## Distorted line profiles
		self.distorted_lines = []
		for ii in range(len(self.spot_rings)):

			## Make the limb-darkened (LD) stellar grid into rings
			LD_ring = self.spot_rings[ii]*self.ring_grid
			self.active_ring = LD_ring.copy()

			## Microturbulence 
			## only micro-turbulence is different
			## WHEN implemented this way
			self.microturbulence(self.xi)

			## Stellar line profile
			self.convolve()

			line = np.add.reduce(self.line_conv,axis=0)
			dline = line_oot - line
			self.distorted_lines.append(dline/self.area)


	def showLines(self,ax=None):

		if ax is None:
			fig = plt.figure()
			ax = fig.add_subplot(111)

		for ii in range(len(self.distorted_lines)):
			ax.plot(self.vel_1d,self.line_profile-self.distorted_lines[ii],color='k',lw=3.0)
			ax.plot(self.vel_1d,self.line_profile-self.distorted_lines[ii],lw=2.0)
		ax.set_xlabel('Velocity (km s$^{-1}$)')
		ax.set_ylabel('CCF')


