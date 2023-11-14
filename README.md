# obscurae
Create stellar line profiles obscured by a transiting planet or/and a spot crossing the disk

### Documentation
[![Documentation Status](https://readthedocs.org/projects/obscurae/badge/?version=latest)](https://obscurae.readthedocs.io/en/latest/?badge=latest)

### Installation
`cd /path/to`

`git clone https://github.com/emilknudstrup/obscurae.git`

`cd /path/to/obscurae`

`python -m pip install .`

or for an editable version

`python -m pip install -e .`

### Example: Distortion of stellar line from transiting planet
```python

# Create a planet object/stellar grid
import matplotlib.pyplot as plt
import obscurae as obsc
rp = 0.1 # planet-to-star radius ratio
planet = obsc.Planet(rp=rp,vsini=10.,zeta=3.0,xi=2.0,cs=[0.8,0.3])
planet.Grid()
planet.Line()
planet.starLight()
planet.starLine()

# Plot the planet's trace across the stellar disk
import numpy as np
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

## calculate true anomaly, f, or rather cos(f) and sin(f)
cosf, sinf = obsc.Dynamics.trueAnomaly(ts,t0,e,p)
## from that calculate x,y positions on stellar disk
xx, yy = obsc.Dynamics.xyPos(cosf,sinf,e,w,a,io,l)

## plot the trace of the planet
planet.transit(trace=(xx,yy))  #path followed by the planet

## plot the occulted stellar disk
planet.Occult(xx,yy)
planet.transit(trace=(xx,yy)) #shows the trace of the planet at each timestamp

# Plot the deformation of the lines
planet.distortLine()
planet.showLines()
plt.show()

```