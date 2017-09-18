import numpy as np
import scipy.interpolate as si
from chisq2sigma import chisq_to_sigma
from argparse import ArgumentParser
import sys
import warnings

class chi2_nh_ih(object):

      def __init__(self, root='gsnpde'):
          self.root = root

      def get_sigma(self):

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

          for i in range(0,len(chi2_osc_nh)):
              chi2_tot_nh.append(chi2_osc_nh[i]+chi2_cos_nh(mnu_nh[i]))

          for i in range(0,len(chi2_osc_ih)):
              chi2_tot_ih.append(chi2_osc_ih[i]+chi2_cos_ih(mnu_ih[i]))

          deltachi2 = min(chi2_tot_nh)-min(chi2_tot_ih)
          deltachi2 = np.abs(deltachi2)
          dof = 1.0

          return chisq_to_sigma(chisq = deltachi2).n_sigma()

if __name__ == '__main__':

   parser = ArgumentParser(prog=sys.argv[0])

   parser.add_argument("-r","--root",
                       dest="root",
                       default='gsnpde',
                       type=str,
                       help="root")

   args = parser.parse_args()

   root = args.root

   chi_squared = chi2_nh_ih(root = root)

   sigma = chi_squared.get_sigma()

   print sigma
