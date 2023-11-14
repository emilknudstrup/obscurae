import matplotlib.pyplot as plt
import numpy as np
from .shadow import *
from .dynamics import *

__all__ = ['test']

def test():
	print('### Testing obscurae ###')
	
	print('\n Creating grid...')
	rp = 0.1 # planet-to-star radius ratio
	planet = Planet(rp=rp,vsini=10.,zeta=3.0,xi=2.0,cs=[0.8,0.3])
	planet.Grid()
	print('\n check')
	print('\n Creating line profile...')
	planet.Line()
	planet.starLight()
	planet.starLine()
	print('\n check')

	## timestamps/observations
	ts = np.array([-0.04,-0.025,-0.01,0.0,0.02,0.05])
	## short-period orbit 
	p = 2.0 #period in days
	t0 = 0 #mid-transit time

	## impact parameter of ~0.4
	io = np.deg2rad(86) #orbital inclination in radians
	a = 6 #semi-major axis in stellar radii

	## circular orbit
	e = 0.0 #eccentricity
	w = np.deg2rad(90) #argument of periapse in radians
	## aligned orbit
	l = np.deg2rad(0) #projected spin-orbit angle/obliquity in radians

	print('\n Calculating true anomaly...')
	## calculate true anomaly, f, or rather cos(f) and sin(f)
	cosf, sinf = Dynamics.trueAnomaly(ts,t0,e,p)
	print('\n check')
	print('\n Calculating x,y positions...')
	## from that calculate x,y positions on stellar disk
	xx, yy = Dynamics.xyPos(cosf,sinf,e,w,a,io,l)
	print('\n check')

	print('\n Calculating transit...')
	## plot the trace of the planet
	planet.transit(trace=(xx,yy))  #path followed by the planet

	## plot the occulted stellar disk
	planet.Occult(xx,yy)
	planet.transit(trace=(xx,yy)) #shows the trace of the planet at each timestamp
	print('\n check')

	print('\n Calculating distorted line profile...')
	# Plot the deformation of the lines
	planet.distortLine()
	planet.showLines()
	print('\n check')
	print('### Testing obscurae complete ###')
	plt.show()

