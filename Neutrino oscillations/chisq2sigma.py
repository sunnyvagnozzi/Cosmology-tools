import scipy.stats as st
from argparse import ArgumentParser
import sys

class chisq_to_sigma(object):

      def __init__(self, chisq=1.0, dof=1.0):

          self.chisq= chisq
          self.dof = dof

      def n_sigma(self):
          pvalue = st.chi2.pdf(self.chisq, self.dof)
          nsigma = st.norm.ppf(1.0-pvalue/2.0)
          return nsigma

if __name__ == '__main__':

   parser = ArgumentParser(prog=sys.argv[0])

   parser.add_argument("-chisq","--chisquare",
                       dest="chisq",
                       default=1.0,
                       type=float,
                       help="Chi-squared value")

   parser.add_argument("-dof","--degreesoffreedom",
                       dest="dof",
                       default=1.0,
                       type=int,
                       help="Number of degrees of freedom")

   args = parser.parse_args()

   chisq = args.chisq
   dof = args.dof

   chi_squared = chisq_to_sigma(chisq = chisq, dof = dof)

   num_sigma = chi_squared.n_sigma()

   print num_sigma
