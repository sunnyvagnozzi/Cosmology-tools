import numpy as np
import pyfits as pf
import os
from astropy.cosmology import FlatLambdaCDM
import astropy.units as un
from argparse import ArgumentParser
import sys

des = 'Code which reads the BOSS galaxy and random catalogs in fits format \
        and returns text files containing RA, DEC, Z, W, N or \
        X, Y, Z, W, N columns (allowing choice of cosmology). \
        Reading is performed by the read() method \
        defined in the Read_BOSS class.'
         
epi = 'For more information or help email sunny.vagnozzi@fysik.su.se.'

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

class Read_BOSS(object):
      
      def __init__(self,
                   gnamenorth = 'galaxy_DR12v5_CMASS_North.fits', 
                   gnamesouth = 'galaxy_DR12v5_CMASS_South.fits', 
                   rnamenorth = 'random0_DR12v5_CMASS_North.fits', 
                   rnamesouth = 'random0_DR12v5_CMASS_South.fits', 
                   zlimits = [0.43, 0.7], 
                   lens = [777202, 40096792], 
                   save_rdzwn = True, 
                   save_xyzwn = True, 
                   use_sys = True, 
                   hubble=70.0, 
                   omegamatter=0.31):

          self.gnamenorth = gnamenorth
          self.gnamesouth = gnamesouth
          self.rnamenorth = rnamenorth
          self.rnamesouth = rnamesouth
          self.zlimits = zlimits
          self.lens = lens
          self.save_rdzwn = save_rdzwn
          self.save_xyzwn = save_xyzwn
          self.use_sys = use_sys
          self.hubble = hubble
          self.omegamatter = omegamatter

      def read(self):

          def radec2cart(rav, decv, chiv):
              xv = np.abs(chiv*np.cos(rav)*np.cos(decv))
              yv = np.abs(chiv*np.sin(rav)*np.cos(decv))
              zv = np.abs(chiv*np.sin(decv))
              return xv, yv, zv

          cosmo = FlatLambdaCDM(H0=self.hubble * un.km / un.s / un.Mpc,
                                Tcmb0=2.725 * un.K, Om0=self.omegamatter)

          zmin, zmax = self.zlimits[0], self.zlimits[1]

          len_gal, len_ran = self.lens[0], self.lens[1]

          print 'Analyzing the galaxies...\n'

          north, south = (pf.open(self.gnamenorth))[1].data, \
                         (pf.open(self.gnamesouth))[1].data

          northcolumns, southcolumns = (pf.open(self.gnamenorth))[1].columns.names, \
                                       (pf.open(self.gnamesouth))[1].columns.names

          print 'Galaxy columns: \n', 'North: ', northcolumns, '\n', \
                'South: ', southcolumns, '\n'

          ran, decn, zn, wfkpn, nn = north['RA'], north['DEC'], north['Z'], \
                                     north['WEIGHT_FKP'], north['NZ']
          wstarn, wseen, wcpn, wnozn = north['WEIGHT_STAR'], north['WEIGHT_SEEING'], \
                                       north['WEIGHT_CP'], north['WEIGHT_NOZ']
          ras, decs, zs, wfkps, ns = south['RA'], south['DEC'], south['Z'], \
                                     south['WEIGHT_FKP'], south['NZ']
          wstars, wsees, wcps, wnozs = south['WEIGHT_STAR'], south['WEIGHT_SEEING'], \
                                       south['WEIGHT_CP'], south['WEIGHT_NOZ']

          indexnorth, indexsouth = (zn >= zmin) & (zn <=zmax), \
                                   (zs >= zmin) & (zs <=zmax)

          ran, decn, zn, wfkpn, nn = ran[indexnorth], decn[indexnorth], zn[indexnorth], \
                                     wfkpn[indexnorth], nn[indexnorth]
          wstarn, wseen, wcpn, wnozn = wstarn[indexnorth], wseen[indexnorth], \
                                       wcpn[indexnorth], wnozn[indexnorth]
          ras, decs, zs, wfkps, ns = ras[indexsouth], decs[indexsouth], zs[indexsouth], \
                                     wfkps[indexsouth], ns[indexsouth]
          wstars, wsees, wcps, wnozs = wstars[indexsouth], wsees[indexsouth], \
                                       wcps[indexsouth], wnozs[indexsouth]

          if self.use_sys:
              print 'Using systematics \n'
              wn = wfkpn*wstarn*wseen*(wcpn+wnozn-1.0)
              ws = wfkps*wstars*wsees*(wcps+wnozs-1.0)
          else:
              print 'Not using systematics, weights are just FKP weights \n'
              wn = wfkpn
              ws = wfkps

          ra_gal, dec_gal, red_gal, w_gal, n_gal = np.concatenate((ran, ras)), \
                                                   np.concatenate((decn, decs)), \
                                                   np.concatenate((zn, zs)), \
                                                   np.concatenate((wn, ws)), \
                                                   np.concatenate((nn, ns))

          print 'Analyzing the randoms, have patience...\n'

          north, south = (pf.open(self.rnamenorth))[1].data, \
                         (pf.open(self.rnamesouth))[1].data

          northcolumns, southcolumns = (pf.open(self.rnamenorth))[1].columns.names, \
                                       (pf.open(self.rnamesouth))[1].columns.names

          print 'Random columns: \n', 'North: ', northcolumns, '\n', \
                'South: ', southcolumns, '\n'

          ran, decn, zn, wfkpn, nn = north['RA'], north['DEC'], north['Z'], \
                                     north['WEIGHT_FKP'], north['NZ']
          ras, decs, zs, wfkps, ns = south['RA'], south['DEC'], south['Z'], \
                                     south['WEIGHT_FKP'], south['NZ']

          indexnorth, indexsouth = (zn >= zmin) & (zn <= zmax), \
                                   (zs >= zmin) & (zs <= zmax)

          ran, decn, zn, wfkpn, nn = ran[indexnorth], decn[indexnorth], zn[indexnorth], \
                                     wfkpn[indexnorth], nn[indexnorth]
          ras, decs, zs, wfkps, ns = ras[indexsouth], decs[indexsouth], zs[indexsouth], \
                                     wfkps[indexsouth], ns[indexsouth]

          ra_ran, dec_ran, red_ran, w_ran, n_ran = np.concatenate((ran, ras)), \
                                                   np.concatenate((decn, decs)), \
                                                   np.concatenate((zn, zs)), \
                                                   np.concatenate((wfkpn, wfkps)), \
                                                   np.concatenate((nn, ns))

          print 'Converting galaxies to desired format...\n'

          print 'Using a fiducial cosmology with Omega_m = ', cosmo.Om(0), \
                ' and h = ', (cosmo.H(0).value)/100.0, '\n'

          ra_gal, dec_gal = np.deg2rad(ra_gal), np.deg2rad(dec_gal)
          chi_gal = cosmo.comoving_distance(red_gal).value
          x_gal, y_gal, z_gal = radec2cart(ra_gal, dec_gal, chi_gal)

          print 'Converting randoms to desired format, have patience...\n'

          print 'Using a fiducial cosmology with Omega_m = ', cosmo.Om(0), \
                ' and h = ', (cosmo.H(0).value)/100.0, '\n'

          ra_ran, dec_ran = np.deg2rad(ra_ran), np.deg2rad(dec_ran)
          chi_ran = cosmo.comoving_distance(red_ran).value
          x_ran, y_ran, z_ran = radec2cart(ra_ran, dec_ran, chi_ran)

          if (len(x_gal)==len_gal and len(y_gal)==len_gal and len(z_gal)==len_gal \
              and len(w_gal)==len_gal and len(n_gal)==len_gal):
             print 'Lengths of galaxies match\n'
          else:
             print 'Lenghts of galaxies do not match :(\n'
             os._exit(1)

          if (len(x_ran)==len_ran and len(y_ran)==len_ran and len(z_ran)==len_ran \
              and len(w_ran)==len_ran and len(n_ran)==len_ran):
             print 'Lengths of randoms match\n'
          else:
             print('Lengths of randoms do not match :(\n')
             os._exit(1)

          if self.save_rdzwn:
             print 'Saving to rdzwn files, have patience...\n'
             np.savetxt('galaxies_rdzwn',np.column_stack((ra_gal, dec_gal, \
                         red_gal, w_gal, n_gal)), fmt='%15.8e', \
                         header='Columns are: RA, DEC, redshift, weights, comoving number density. Using a fiducial cosmology with Omegamatter = '+str(self.omegamatter)+' and H0 = '+str(self.hubble)+' km/s/Mpc')
             np.savetxt('randoms_rdzwn', np.column_stack((ra_ran, dec_ran, \
                         red_ran, w_ran, n_ran)), fmt='%15.8e',
                         header='Columns are: RA, DEC, redshift, weights, comoving number density. Using a fiducial cosmology with Omegamatter = '+str(self.omegamatter)+' and H0 = '+str(self.hubble)+' km/s/Mpc')
             print 'Saved to the following files: \n', \
                   'galaxies_rdzwn and randoms_rdzwn \n'

          if self.save_xyzwn:
             print 'Saving to xyzwn files, have patience...\n'
             np.savetxt('galaxies_xyzwn', np.column_stack((x_gal, y_gal, \
                        z_gal, w_gal, n_gal)), fmt='%15.8e',
                        header='Columns are: Cartesian X, Y, Z, weights, comoving number density. Using a fiducial cosmology with Omegamatter = '+str(self.omegamatter)+' and H0 = '+str(self.hubble)+' km/s/Mpc')
             np.savetxt('randoms_xyzwn', np.column_stack((x_ran, y_ran, \
                        z_ran, w_ran, n_ran)), fmt='%15.8e',
                        header='Columns are: Cartesian X, Y, Z, weights, comoving number density. Using a fiducial cosmology with Omegamatter = '+str(self.omegamatter)+' and H0 = '+str(self.hubble)+' km/s/Mpc')
             print 'Saved to the following files: \n', \
                   'galaxies_xyzwn and randoms_rdzwn \n'

if __name__ == '__main__':
   
   parser = ArgumentParser(prog=sys.argv[0], \
                           description=des, \
                           epilog=epi)

   parser.add_argument("-gn","--galnorth",
                        dest="gnamenorth",
                        default='galaxy_DR12v5_CMASS_North.fits',
                        type=str,
                        help="name of Northern galaxy catalogue (with extension)")

   parser.add_argument("-gs","--galsouth",
                        dest="gnamesouth",
                        default="galaxy_DR12v5_CMASS_South.fits",
                        type=str,
                        help="name of Southern galaxy catalogue (with extension)")

   parser.add_argument("-rn","--rannorth",
                        dest="rnamenorth",
                        default="random0_DR12v5_CMASS_North.fits",
                        type=str,
                        help="name of Northern random catalogue (with extension)")

   parser.add_argument("-rs","--ransouth",
                        dest="rnamesouth",
                        default="random0_DR12v5_CMASS_South.fits",
                        type=str,
                        help="name of Southern random catalogue (with extension)")

   parser.add_argument("-z","--redshift",
                        dest="zlimits",
                        default=[0.43, 0.7],
                        type = list,
                        help="upper and lower limit for redshift cut \
                              (given as a list, default is [0.43, 0.7])")

   parser.add_argument("-l","--length",
                        dest="lens",
                        default=[777202, 40096792],
                        type=list,
                        help="number of galaxies and randoms after cut \
                              (given as a list, default is [777202, 40096792])")

   parser.add_argument("-r","--rdzwn",
                        dest="save_rdzwn",
                        default=True,
                        type=str2bool,
                        help="whether or not to save RA, DEC, Z, W, N files \
                             (answer Yes, True', T, Y, 1 or No, False, F, N, 1)")

   parser.add_argument("-x","--xyzwn",
                        dest="save_xyzwn",
                        default=True,
                        type=str2bool,
                        help="whether or not to save X, Y, Z, W, N files \
                             (answer Yes, True', T, Y, 1 or No, False, F, N, 1)")

   parser.add_argument("-s","--sys",
                        dest="use_sys",
                        default=True,
                        type=str2bool,
                        help="whether or not to use systematics for the galaxy catalogue \
                             (answer Yes, True', T, Y, 1 or No, False, F, N, 1)")

   parser.add_argument("-hub","--hubble",
                        dest="hubble",
                        default=70.0,
                        type=float,
                        help="Hubble parameter, H0, in units of km/s/Mpc")

   parser.add_argument("-m","--matter",
                        dest="omegamatter",
                        default=0.31,
                        type=float,
                        help="Omega_matter")

   args = parser.parse_args()

   gnamenorth = args.gnamenorth
   gnamesouth = args.gnamesouth
   rnamenorth = args.rnamenorth
   rnamesouth = args.rnamesouth
   zlimits = args.zlimits
   lens = args.lens
   save_rdzwn = args.save_rdzwn
   save_xyzwn = args.save_xyzwn
   use_sys = args.use_sys
   hubble = args.hubble
   omegamatter = args.omegamatter

   data = Read_BOSS(gnamenorth = gnamenorth,
                    gnamesouth = gnamesouth,
                    rnamenorth = rnamenorth,
                    rnamesouth = rnamesouth,
                    zlimits = zlimits,
                    lens = lens,
                    save_rdzwn = save_rdzwn,
                    save_xyzwn = save_xyzwn,
                    use_sys = use_sys,
                    hubble = hubble,
                    omegamatter = omegamatter)

   data.read()
