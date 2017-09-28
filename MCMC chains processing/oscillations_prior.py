import numpy as np
from argparse import ArgumentParser
import sys
import os

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False

class oscillations_prior(object):
      
      def __init__(self,
                   root = 'pgg_clkg_2params_rsd_bk_mnu',
                   num = 6,
                   param_index = 5,
                   param_limit = 0.06,
                   greater = True):

          self.root = root
          self.num = num
          self.param_index = param_index
          self.param_limit = param_limit
          self.greater = greater

      def impose_prior(self):

          for i in range(1,self.num+1):
              num_chain = str(i)
              chain = np.loadtxt(self.root+'_'+num_chain+'.txt')
              param = chain[:,self.param_index+1]

              if self.greater:
                 index = param >= self.param_limit

                 if i==1:
                    print "Processing chains using condition param >=", self.param_limit
              else:
                 index = param <= self.param_limit

                 if i==1:
                    print "Processing chains using condition param <=", self.param_limit

              chain_new = chain[index,:]
              np.savetxt(self.root+'_oscillations_'+num_chain+'.txt', chain_new, fmt='%15.8e')
              print "Processed chain", num_chain

          command = 'cp '+self.root+'.paramnames '+self.root+'_oscillations.paramnames'
          os.system(command)
          print "Copied paramnames"

if __name__ == '__main__':
   
   parser = ArgumentParser(prog=sys.argv[0],
                           add_help = True)

   parser.add_argument("-r","--root",
                       dest="root",
                       default='pgg_clkg_2params_rsd_bk_mnu',
                       type=str,
                       help="Root of chains ([root]_[chain_num].txt). Default: 'pgg_clkg_2params_rsd_bk_mnu'")

   parser.add_argument("-n","--num",
                       dest="num",
                       default=6,
                       type=int,
                       help="Number of chains to process. Default: 6")

   parser.add_argument("-pi","--parami",
                       dest="param_index",
                       default=5,
                       type=int,
                       help="Index of parameter to be processed. Default: 5 (mnu)")

   parser.add_argument("-pl","--paraml",
                       dest="param_limit",
                       default=0.06,
                       type=float,
                       help="Limit of new prior on paramter. Default: 0.06")

   parser.add_argument("-g","--greater",
                       dest="greater",
                       default=True,
                       type=str2bool,
                       help="Whether to process the chains using param >= param_limit or param <= param_limit\
                             (answer Yes, True, T, Y, 1 or No, False, F, N, 1, default is False). Default: True")


   args = parser.parse_args()

   root = args.root
   num = args.num
   param_index = args.param_index
   param_limit = args.param_limit
   greater = args.greater

   old_chains = oscillations_prior(root = root,
                                   num = num,
                                   param_index = param_index,
                                   param_limit = param_limit,
                                   greater = greater)

   old_chains.impose_prior()
