import os

command1 = 'python DR11reader.py'
command2 = 'python DR11formatter.py'
command3 = 'cc fkpDR11.cpp -o fkpDR11 -lm -lfftw3_threads -lfftw3 -lpthread -fopenmp'
command4 = './fkpDR11 -gal galaxies_DR11 -ran randoms_DR11 -galext dat -ranext dat -out psDR11 -outext dat'

for number in range(1,5):
    if number == 1:
       print 'Reading DR11 galaxies and randoms, printing columns, implementing CMASS cut (0.43 < z < 0.7) and combining North and South catalogs \n'
       os.system(command1)
       print 'Reading DR11 successful and combining North and South catalogs successful. Output is in galaxies_DR11_CMASS_NS.dat and randoms_DR11_CMASS_NS.dat \n'
    elif number == 2:
       print 'Putting galaxies and randoms in correct format: (x,y,z,w,n). Printing comoving distance every 100000 objects \n'
       os.system(command2)
       print 'Formatting executed successfully. Output is in galaxies_DR11.dat and randoms_DR11.dat \n'
    elif number == 3:
       print 'Compiling power spectrum code \n'
       os.system(command3)
       print 'Power spectrum code compiled successfully. Executable is fkpDR11 \n'
    elif number == 4:
       print 'Executing power spectrum code \n'
       os.system(command4)
       print 'Power spectrum code executed successfully. Output is in psDR11.dat \n'

print 'DR11 power spectrum pipeline executed successfully \n'
