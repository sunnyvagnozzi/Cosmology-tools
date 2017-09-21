import numpy as np
from argparse import ArgumentParser
import sys

class omega(object):

      def __init__(self,
                   root = 'gsnpde',
                   num = 8):
          
          self.root = root
          self.num = num

      def process_chains(self):

          for i in range(1,self.num+1):
              num_chain = str(i)
              chain = np.loadtxt(self.root+'_'+num_chain+'.txt')
              omegabh2, omegach2, H0 = chain[:,2], chain[:,3], chain[:,40]
              h = H0/100.0
              omega = (omegabh2 + omegach2)/(h*h)
              chain_new = np.column_stack((chain, omega))
              np.savetxt(self.root+'_omega_'+str(i)+'.txt', chain_new, fmt='%15.8e')
              print "Processed chain ", num_chain
              del chain
              del omegabh2
              del omegach2
              del H0
              del h
              del omega
              del chain_new

if __name__ == '__main__':

   parser = ArgumentParser(prog=sys.argv[0],
                           add_help = True)

   parser.add_argument("-r","--root",
                       dest="root",
                       default='gsnpde',
                       type=str,
                       help="Root of chains ([root]_[chain_num].txt). Default: 'gsnpde'")

   parser.add_argument("-n","--num",
                       dest="num",
                       default=8,
                       type=int,
                       help="Number of chains to process. Default: 8")

   args = parser.parse_args()

   root = args.root
   num = args.num

   omega(root = root, num = num).process_chains()
