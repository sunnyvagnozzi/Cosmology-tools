import numpy as np
import pyfits as pf
import os

rag = []
decg = []
redg = []
wg = []
ng = []

rar = []
decr = []
redr = []
wr = []
nr = []

for i in range(0,4):
    datatype = i

    if datatype == 0:
       fname = '../Data/galaxy_DR11v1_CMASS_North'
       datafile = pf.open(fname + '.fits')
    if datatype == 1:
       fname = '../Data/galaxy_DR11v1_CMASS_South'
       datafile = pf.open(fname + '.fits')
    if datatype == 2:
       fname = '../Data/random0_DR11v1_CMASS_North'
       datafile = pf.open(fname + '.fits')
    if datatype == 3:
       fname = '../Data/random0_DR11v1_CMASS_South'
       datafile = pf.open(fname + '.fits')

    data = datafile[1].data

    datacols = datafile[1].columns
    print datacols.names

    ra = data['RA']
    dec = data['DEC']
    red = data['Z']
    w = data['WEIGHT_FKP']
    n = data['NZ']

    count = 0

    for j in range(0, len(red)):
        if red[j] >= 0.43 and red[j] <= 0.7:
           if i == 0:
              rag.append(ra[j])
              decg.append(dec[j])
              redg.append(red[j])
              wg.append(w[j])
              ng.append(n[j])
              count = count + 1
           if i == 1:
              rag.append(ra[j])
              decg.append(dec[j])
              redg.append(red[j])
              wg.append(w[j])
              ng.append(n[j])
              count = count + 1
           if i == 2:
              rar.append(ra[j])
              decr.append(dec[j])
              redr.append(red[j])
              wr.append(w[j])
              nr.append(n[j])
              count = count + 1
           if i == 3:
              rar.append(ra[j])
              decr.append(dec[j])
              redr.append(red[j])
              wr.append(w[j])
              nr.append(n[j])
              count = count + 1

    if i == 0:
        print 'The number of North galaxies surviving the cut is: ', count
    elif i == 1:
        print 'The number of South galaxies surviving the cut is: ', count
    elif i == 2:
        print 'The number of North randoms surviving the cut is: ', count
    elif i == 3:
        print 'The number of South randoms surviving the cut is: ', count

    galaxy = np.column_stack((rag, decg, redg, wg, ng))
    random = np.column_stack((rar, decr, redr, wr, nr))

    np.savetxt('../Outputs/galaxies_DR11_CMASS_NS.dat', galaxy, fmt = '%15.8e')
    np.savetxt('../Outputs/randoms_DR11_CMASS_NS.dat', random, fmt = '%15.8e')
