import numpy as np
import os
from numpy import *

omegam = 0.307115
omegal = 1.0 - omegam

stepsize = 0.001
n = 10
ch0 = 3000.0

def chi(z):
    DCMR = 0.0

    az = 1.0/(1.0 + 1.0*z)

    for i in range(n):
        a = az + (1.0 - az)*(i + 0.5)/n
        adot = np.sqrt((omegam/a) + (omegal*a*a))
        DCMR = DCMR + 1.0/(a*adot)

    DCMR = (1.0 - az)*DCMR/n
    DCMR = DCMR*ch0

    return DCMR;

ragalaxy = []
decgalaxy = []
redgalaxy = []
wgalaxy = []
ngalaxy = []

rarandom = []
decrandom = []
redrandom = []
wrandom = []
nrandom = []

gal = open('../Outputs/galaxies_DR11_CMASS_NS.dat')
ran = open('../Outputs/randoms_DR11_CMASS_NS.dat')

print 'Using a fiducial cosmology with Omega_m = ', omegam, ' and Omega_L = ', 1.0 - omegam

for i in gal:
    g = map(float, i.split())
    ragalaxy.append(g[0])
    decgalaxy.append(g[1])
    redgalaxy.append(g[2])
    wgalaxy.append(g[3])
    ngalaxy.append(g[4])

print 'Number of galaxies: ', len(ragalaxy)

alphagalaxy = np.deg2rad(ragalaxy)
deltagalaxy = np.deg2rad(decgalaxy)

chig = np.zeros(len(redgalaxy))
for i1, i2 in enumerate(redgalaxy):
    chig[i1] = chi(redgalaxy[i1])
    if i1%100000 == 0:
       print i1, chig[i1]

xgalaxy = np.absolute(chig*np.cos(alphagalaxy)*np.cos(deltagalaxy))
ygalaxy = np.absolute(chig*np.sin(alphagalaxy)*np.cos(deltagalaxy))
zgalaxy = np.absolute(chig*np.sin(deltagalaxy))

galaxy = np.column_stack((xgalaxy, ygalaxy, zgalaxy, wgalaxy, ngalaxy))
np.savetxt('../Outputs/galaxies_DR11.dat', galaxy, fmt = '%15.8e')

for j in ran:
    r = map(float, j.split())
    rarandom.append(r[0])
    decrandom.append(r[1])
    redrandom.append(r[2])
    wrandom.append(r[3])
    nrandom.append(r[4])

print 'Number of randoms: ', len(rarandom)

alpharandom = np.deg2rad(rarandom)
deltarandom = np.deg2rad(decrandom)

chir = np.zeros(len(redrandom))
for j1, j2 in enumerate(redrandom):
    chir[j1] = chi(redrandom[j1])
    if j1%100000 == 0:
       print j1, chir[j1]

xrandom = np.absolute(chir*np.cos(alpharandom)*np.cos(deltarandom))
yrandom = np.absolute(chir*np.sin(alpharandom)*np.cos(deltarandom))
zrandom = np.absolute(chir*np.sin(deltarandom))

random = np.column_stack((xrandom, yrandom, zrandom, wrandom, nrandom))
np.savetxt('../Outputs/randoms_DR11.dat', random, fmt = '%15.8e')
