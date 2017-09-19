#!usr/bin/env python
"""
Version : 0.1.0
Date : 19th September 2017
Authors : Sunny Vagnozzi
Email : sunny.vagnozzi@fysik.su.se
Affiliation : The Oskar Klein Centre for Cosmoparticle Physics - Stockholm University
License : MIT
Status : Under Development
Description :
Code which combines posteriors on the sum of the three active neutrino masses obtained from cosmological data
with the Chi^2 profiles for the solar and atmospheric mass splittings from NuFIT 3.0 (http://www.nu-fit.org/?q=node/139)
to output the significance (provided as # of sigmas) at which the normal neutrino mass ordering is preferred
over the inverted neutrino mass ordering.
The code provides the chi2_neutrino_hierarchy class, wherein the num_sigma() method calculates the significance
"""

import numpy as np
import scipy.interpolate as si
from chisq2sigma import chisq_to_sigma
from argparse import ArgumentParser
import sys
import warnings

__author__ = "Sunny Vagnozzi"
__email__ = "sunny.vagnozzi@fysik.su.se"
__license__ = "MIT"
__version__ = "0.1.0"
__status__ = "Development"

des = '**\nCode which combines cosmology and neutrino oscillation data from NuFIT 3.0 \
       to determine the # of sigmas at which the normal hierarchy is preferred \
       over the inverted hierarchy.\n**\n'

epi = '**\nDeveloped by Sunny Vagnozzi, 2017. \
       For more information, help, or for reporting bugs please contact sunny.vagnozzi@fysik.su.se. \
       For help with the usage of the code, type in the command line: python chi2_neutrino_hierarchy.py -h \
       Known bugs: when running with the flag "-d 2" (i.e. when only cosmological data is usecd) \
       the code will return nan sigma since cosmological data is not sensitive to the neutrino mass ordering \
       and hence Delta Chi^2 (NH-IH) = 0. \
       Will be improved in a future release through the treatment of Hannestad & Schwetz 2016.\n**\n'

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False

class chi2_neutrino_hierarchy(object):

      def __init__(self, root = 'gsnpde',
                   save_profiles = True,
                   dataset = 0,
                   verbose = 0):

          self.root = root
          self.save_profiles = save_profiles
          self.dataset = dataset
          self.verbose = verbose

      def num_sigma(self):

          mnu_nh = []
          mnu_ih = []

          mnu_cos_nh = []
          mnu_cos_ih = []

          chi2_osc_nh = []
          chi2_osc_ih = []

          chi2_cos_nh = []
          chi2_cos_ih = []

          chi2_tot_nh = []
          chi2_tot_ih = []

          deltam2_21_nh, chi2_21_nh = np.loadtxt('deltam2_21_nh.dat', unpack=True)
          deltam2_31_nh, chi2_31_nh = np.loadtxt('deltam2_31_nh.dat', unpack=True)
          deltam2_21_ih, chi2_21_ih = np.loadtxt('deltam2_21_ih.dat', unpack=True)
          deltam2_32_ih, chi2_32_ih = np.loadtxt('deltam2_32_ih.dat', unpack=True)

          for i in range(0,len(deltam2_31_nh)):
              for j in range(0,len(deltam2_21_nh)):
                  mnu_nh.append(np.sqrt(deltam2_31_nh[i]*0.001) \
                               +np.sqrt(10.0**(deltam2_21_nh[j])))
                  chi2_osc_nh.append(chi2_31_nh[i]+chi2_21_nh[j])

          for i in range(0,len(deltam2_32_ih)):
              for j in range(0,len(deltam2_21_ih)):
                  mnu_try = np.sqrt(-deltam2_32_ih[i]*0.001) \
                           +np.sqrt((-10.0**(deltam2_21_ih[j]))-deltam2_32_ih[i]*0.001)
                  warnings.filterwarnings("ignore")
                  if np.isnan(mnu_try)==0:
                     mnu_ih.append(mnu_try)
                     chi2_osc_ih.append(2.0*chi2_32_ih[i]+chi2_21_ih[j])

          mnu_nh = np.array(mnu_nh)
          chi2_osc_nh = np.array(chi2_osc_nh)

          mnu_ih = np.array(mnu_ih)
          chi2_osc_ih = np.array(chi2_osc_ih)

          mnu, pos = np.loadtxt(self.root+'_p_mnu.dat', unpack=True)

          chi2_cos = -2.0*np.log(pos)

          chi2_cos_nh = si.interp1d(mnu,chi2_cos)
          chi2_cos_ih = si.interp1d(mnu,chi2_cos)

          if self.dataset == 0 and self.verbose > 0:
             print 'Using data from both neutrino oscillations and cosmology\n'
          elif self.dataset == 1 and self.verbose > 0:
             print 'Using only data from neutrino oscillations\n'
          elif self.verbose > 0:
             print 'Using only data from cosmology\n'

          for i in range(0,len(chi2_osc_nh)):
              if self.dataset == 0:
                 chi2_tot_nh.append(chi2_osc_nh[i]+chi2_cos_nh(mnu_nh[i]))
              elif self.dataset == 1:
                 chi2_tot_nh.append(chi2_osc_nh[i])
              else:
                 chi2_tot_nh.append(chi2_cos_nh(mnu_nh[i]))

          for i in range(0,len(chi2_osc_ih)):
              if self.dataset == 0:
                 chi2_tot_ih.append(chi2_osc_ih[i]+chi2_cos_ih(mnu_ih[i]))
              elif self.dataset == 1:
                 chi2_tot_ih.append(chi2_osc_ih[i])
              else:
                 chi2_tot_ih.append(chi2_cos_nh(mnu_ih[i]))

          if self.save_profiles:
             if self.verbose > 0:
                print 'Saving Chi^2 profiles...\n'
                chi2_nh, chi2_ih = np.column_stack((mnu_nh, chi2_tot_nh)), \
                                   np.column_stack((mnu_ih, chi2_tot_ih))

                np.savetxt('profile_chi2_nh.dat', chi2_nh, fmt='%15.8e')
                np.savetxt('profile_chi2_ih.dat', chi2_ih, fmt='%15.8e')
                print 'Chi^2 profiles saved to '
                print 'profile_chi2_nh.dat (for normal hierarchy) '
                print 'and '
                print 'profile_chi2_ih.dat (for inverted hierarchy)\n'
             

          deltachi2 = min(chi2_tot_nh)-min(chi2_tot_ih)
          deltachi2 = np.abs(deltachi2)
          dof = 1.0

          return round(chisq_to_sigma(chisq = deltachi2).n_sigma(),1)

if __name__ == '__main__':

   parser = ArgumentParser(prog=sys.argv[0],
                           add_help = True,
                           description = des,
                           epilog = epi)

   parser.add_argument("-r","--root",
                       dest="root",
                       default='gsnpde',
                       type=str,
                       help="Root of the M_nu posterior file ([root]_p_mnu.dat),\
                             includes the path here if required (default: 'gsnpde').")

   parser.add_argument("-s","--save",
                        dest="save_profiles",
                        default=False,
                        type=str2bool,
                        help="Whether or not to save Chi^2 profiles against M_nu \
                             (answer Yes, True, T, Y, 1 or No, False, F, N, 1, default is False)")

   parser.add_argument("-d","--dataset",
                        dest="dataset",
                        default=0,
                        type=int,
                        help="Which datasets to use to calculate the Chi^2\
                             0: oscillations and cosmology (default)\
                             1: only oscillations\
                             2: only cosmology")

   parser.add_argument("-v","--verbose",
                        dest="verbose",
                        default=0,
                        type=int,
                        help="Controls chattiness of code:\
                              0 (default) for non chatty, >0 for chatty.")

   args = parser.parse_args()

   root = args.root
   save_profiles = args.save_profiles
   dataset = args.dataset
   verbose = args.verbose

   chi2_profile = chi2_neutrino_hierarchy(root = root,
                            save_profiles = save_profiles,
                            dataset = dataset,
                            verbose = verbose)

   sigma = chi2_profile.num_sigma()

   print '\nNormal hierarchy preferred over inverted hierarchy at ', sigma,u'\u03C3.\n'
