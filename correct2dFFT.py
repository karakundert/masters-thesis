#! /usr/bin/env python
from numpy import *
from scipy import *
from scipy import fftpack
from scipy import signal
from numpy import fft

import matplotlib.pylab as pl

########## Unit conversion
##  1.0 / ( cell size in radians ) = ( u_max in units of wavelength )
##  1.0 / ( field of view in radians ) = ( cell size in uv-domain in units of wavelength )

cell = 0.1 # arcmin
fov = 30.0 # arcmin
freq = 1.5e9
c = 3e8
wvlen = c / freq # observed wavelength
d = 100 # dish width
spat_lam = d / wvlen


# Image pixels
xvals = arange(-1*fov/2.0,fov/2.0,cell)
yvals = arange(-1*fov/2.0,fov/2.0,cell)
Nxy = xvals.shape[0]

# UV cell size and extent - in units of wavelengths
uvcell = 1 / (fov/60.0 * pi/180.0)
uvmax = Nxy/2 * uvcell

uvals = arange(-1*uvmax, uvmax, uvcell)
vvals = arange(-1*uvmax, uvmax, uvcell)
Nuv = uvals.shape[0]

print 'cell (arcmin): ', cell
print 'fov (arcmin): ', fov
print 'uvcell (lambda): ',uvcell
print 'uvmax (lambda): ', uvmax

print Nxy, Nuv  # To make sure they're equal !

###################

# Aperture Function
aper = zeros((Nuv,Nuv))
d_uv = spat_lam / uvcell
uu, vv = mgrid[:Nuv, :Nuv]
circle = (uu - (Nuv/2.0)) ** 2 + (vv - (Nuv/2.0)) ** 2
disk = circle < (4*d_uv/2.0)**2
aper[disk] = 1.0
for u in xrange(Nuv):
    for v in xrange(Nuv):
        if disk[u][v] == True:
            aper[u][v] = 1.0

bessel = (Nxy**2) * fftpack.ifft2(fftpack.fftshift(aper**2))
bessel = fftpack.ifftshift(bessel)

pl.figure(1)
pl.clf()
pl.subplot(121)
pl.imshow(abs(aper**2)[100:200,100:200],extent=[-50,50,-50,50])
pl.subplot(122)
pl.imshow(abs(bessel)[100:200,100:200],extent=[-50,50,-50,50])
pl.savefig('top-hat-bessel.png')

pl.show()
