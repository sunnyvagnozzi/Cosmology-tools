import numpy
from numpy import *
tauold = 0.058
sigmaold = 0.012
taunew = 0.055
sigmanew = 0.009

print "Enter the chain root"
root = raw_input()

for i in range (1,7):
    numchain = str(i)
    mychain = loadtxt(root+'_'+numchain+'.txt')
    mychain[:,0] = mychain[:,0]*(numpy.exp(-((mychain[:,5]-taunew)/sigmanew)**(2.0))/numpy.exp(-((mychain[:,5]-tauold)/sigmaold)**(2.0)))
    numpy.savetxt(root+'_'+numchain+'.txt',mychain,fmt='%15.8e')
    del mychain
