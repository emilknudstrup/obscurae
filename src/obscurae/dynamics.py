
# =============================================================================
# Packages
# =============================================================================

import numpy as np
import scipy.signal as ss

# =============================================================================
# Orbital dynamics
# =============================================================================

class Dynamics(object):
	'''Orbital dynamics.

	Class to calculate the orbital dynamics of the planet,
	i.e., the position of the planet on the stellar disk.

	'''

	@staticmethod
	def keplersEq(mean_anomaly, ecc, tolerance=1.e-5, iter_max=300):
		'''Solves Kepler's equation.

		Function that solves Kepler's equation:
		.. :math:`M = E - \sin(E)`,
		
		where :math:`M` is the mean anomaly and :math:`E` the eccentric anomaly.

		This is done following the Newton-Raphson method as described in :cite:t:`Murray2010`.

		:param mean_anomaly: The mean anomaly.
		:type mean_anomaly: array
		
		:param ecc: Eccentricity.
		:type ecc: float

		:param tolerance: The tolerance for convergene. Defaults to 1.e-5.
		:type tolerance: float, optional

		:param iter_max: Maximum number of iterations. Defaults to 300.
		:type iter_max: int, optional

		:return: The new eccentric anomaly.
		:rtype: array 

		.. note::
			Test if scipy.optimize.newton (https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html) is faster.
		
		'''
		## Circular orbit
		if ecc == 0: return mean_anomaly 

		new_ecc_anomaly = mean_anomaly
		converged = False

		for ii in range(iter_max):
			old_ecc_anomaly = new_ecc_anomaly

			new_ecc_anomaly = old_ecc_anomaly - (old_ecc_anomaly - ecc*np.sin(old_ecc_anomaly) - mean_anomaly)/(1.0 - ecc*np.cos(old_ecc_anomaly))

			if np.max(np.abs(new_ecc_anomaly - old_ecc_anomaly)/old_ecc_anomaly) < tolerance:
				converged = True
				break

		if not converged:
			print('Calculation of the eccentric anomaly did not converge!')

		return new_ecc_anomaly

	@staticmethod
	def trueAnomaly(time, Tw, ecc, per):
		'''Function that returns the true anomaly.

		The approach follows :cite:t:`Murray2010`.
		

		:param time: Times of observations.
		:type time: array

		:param Tw: Time of periastron.
		:type Tw: float

		:param ecc: Eccentricity.
		:type ecc: float

		:param per: Orbital period.
		:type per: float

		:param ww: Argument of periastron in radians.
		:type ww: float


		:return: cosine, sine of the true anomaly.
		:rtype: (array, array)

		'''
		
		n = 2.0*np.pi/per
		mean_anomaly = n*(time-Tw)
		ecc_anomaly = Dynamics.keplersEq(mean_anomaly,ecc)

		cos_E = np.cos(ecc_anomaly)
		sin_E = np.sin(ecc_anomaly)

		## Cosine and sine of the true anomaly
		cos_f = (cos_E - ecc)/(1.0 - ecc*cos_E)
		sin_f = (np.sqrt(1 - ecc**2)*sin_E)/(1.0 - ecc*cos_E)

		return cos_f, sin_f		

	@staticmethod
	def xyPos(cos_f,sin_f,ecc,ww,ar,inc,lam):
		'''Position of planet on stellar disk.

		Function to calculate the position on the stellar disk.
		Stellar disk goes from 0 to 1 in x and y.

		:param cos_f: cosine of the true anomaly
		:type cos_f: array
		:param sin_f: sine of the true anomaly
		:type sin_f: array            
		:param ecc: Eccentricity.
		:type ecc: float
		:param ww: Argument of periastron in radians.
		:type ww: float
		:param ar: Semi-major axis in stellar radii.
		:type ar: float
		:param inc: Inclination in radians.
		:type inc: float
		:param lam: Projected obliquity in radians.
		:type lam: float

		:return: x,y position of planet on stellar disk.
		:rtype: (array, array)
		

		'''
		r = ar*(1.0 - ecc**2)/(1.0 + ecc*cos_f)
		f = np.arctan2(sin_f,cos_f)
		
		## x and y are lists of the positions of the transitting planet on the stellar disk 
		## normalized to stellar radius (using a/Rs), corresponding to each RV-point
		x_old = -1*r*np.cos(ww + f)
		y_old = -1*r*np.sin(ww + f)*np.cos(inc)

		## Rotate our coordinate system, such that the projected obliquity becomes the new y-axis
		x = x_old*np.cos(lam) - y_old*np.sin(lam)
		y = x_old*np.sin(lam) + y_old*np.cos(lam)
		return x, y



# =============================================================================
# Spot position
# =============================================================================

class SpotOn(object):
	'''Spot position.

	Class to calculate the position of a spot on the stellar disk.

	'''

	@staticmethod
	def spotPos(time, theta, phi, per, t_ref=0.0):
		'''Position of spot on stellar surface.

		Function to calculate the position of a spot at a given 
		(:math:`\phi`,:math:`\theta`) on the stellar surface,
		where :math:`\phi` is the colatitude 
		(:math:`0^\circ` at north pole and :math:`180^\circ` at the south pole) 
		and :math:`\theta` the longitude.

		:param time: Times of observations in days.
		:type time: array

		:param phi: Latitude of spot in radians.
		:type phi: float

		:param theta: Longitude of spot in radians.
		:type theta: float

		:param per: Period of spot crossing/rotation period in days.
		:type per: float

		:param t_ref: Reference time in days. Defaults to 0.0.
		:type t_ref: float, optional

		:return: x, y position of spot on stellar surface.
		:rtype: (array, array)
		
		'''
		rot_phase = 2*np.pi*(time-t_ref)/per

		
		## Spot position as a spherical cap on the stellar surface
		## as a function of time
		## note that the z-axis is along the line of sight
		x = np.sin(phi)*np.sin(theta + rot_phase)
		y = np.cos(phi)*np.ones_like(time)


		return x, y


