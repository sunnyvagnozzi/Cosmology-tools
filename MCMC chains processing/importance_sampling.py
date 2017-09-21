import numpy as np
from argparse import ArgumentParser
import sys

class importance_sampling(object):
      
      def __init__(self,
                   root = 'base',
                   num = 6,
                   param_index = 4,
                   param_old = 0.07,
                   sigma_old = 0.02,
                   param_new = 0.05,
                   sigma_new = 0.01):

          self.root = root
          self.num = num
          self.param_index = param_index
          self.param_old = param_old
          self.sigma_old = sigma_old
          self.param_new = param_new
          self.sigma_new = sigma_new

      def process_chains(self):

          for i in range(1,self.num+1):
              num_chain = str(i)
              chain = np.loadtxt(self.root+'_'+num_chain+'.txt')
              counts, param = chain[:,0], chain[:,self.param_index+1]
              counts = counts*\
                      (np.exp(-((param - self.param_new)/self.sigma_new)**(2.0))/\
                       np.exp(-((param - self.param_old)/self.sigma_old)**(2.0)))
              chain[:,0] = counts
              np.savetxt(self.root+'_importance_'+num_chain+'.txt',chain,fmt='%15.8e')
              print "Processed chain ", num_chain
              del chain

if __name__ == '__main__':

   parser = ArgumentParser(prog=sys.argv[0],
                           add_help = True)

   parser.add_argument("-r","--root",
                       dest="root",
                       default='base',
                       type=str,
                       help="Root of chains ([root]_[chain_num].txt). Default: 'base'")

   parser.add_argument("-n","--num",
                       dest="num",
                       default=6,
                       type=int,
                       help="Number of chains to process. Default: 6")

   parser.add_argument("-pi","--parami",
                       dest="param_index",
                       default=4,
                       type=int,
                       help="Index of parameter to process\
                             (\Omega_b*h^2=1, \Omega_c*h^2=2, etc.).\
                             Default: 4 (tau)")

   parser.add_argument("-po","--paramo",
                       dest="param_old",
                       default=0.07,
                       type=float,
                       help="Old value of parameter. Default: 0.07")

   parser.add_argument("-so","--sigmao",
                       dest="sigma_old",
                       default=0.02,
                       type=float,
                       help="Old value of error on parameter. Default: 0.02")

   parser.add_argument("-pn","--paramn",
                       dest="param_new",
                       default=0.05,
                       type=float,
                       help="New value of parameter. Default: 0.05")

   parser.add_argument("-sn","--sigman",
                       dest="sigma_new",
                       default=0.01,
                       type=float,
                       help="New value of error on parameter. Default: 0.01")

   args = parser.parse_args()

   root = args.root
   num = args.num
   param_index = args.param_index
   param_old = args.param_old
   sigma_old = args.sigma_old
   param_new = args.param_new
   sigma_new = args.sigma_new

   post_processing = importance_sampling(root = root,
                                         num = num,
                                         param_index = param_index,
                                         param_old = param_old,
                                         sigma_old = sigma_old,
                                         param_new = param_new,
                                         sigma_new = sigma_new)

   post_processing.process_chains()
