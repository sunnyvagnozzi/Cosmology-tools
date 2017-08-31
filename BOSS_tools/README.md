# BOSS_tools

A set of tools which can be useful to analyze the SDSS-BOSS galaxy (and correspective random) catalogs, especially for first-time users.

## ‬★Read_BOSS

Code which reads the BOSS galaxy and random (North and South) fits files, relying on the read() method defined within the Read_BOSS class. The following 4 outputs are generated:

-*galaxies_rdzwn*: text file with 5 columns containing RA, DEC, redshift, weights and number density of galaxies which pass the redshift cut

-*galaxies_xyzwn*: text file with 5 columns containing X, Y, Z cartesian coordinates, weights and number density of galaxies which pass the redshift cut

-*randoms_rdzwn*: text file with 5 columns containing RA, DEC, redshift, weights and number density of randoms which pass the redshift cut

-*randoms_xyzwn*: text file with 5 columns containing X, Y, Z cartesian coordinates, weights and number density of randoms which pass the redshift cut

### Customizing the code

There are a number of variables

### Usage

### To run Read_BOSS from an ipython terminal or notebook

To run Read_BOSS from an ipython terminal or notebook (in the same folder as Read_BOSS.py) simply type:

    >> from Read_BOSS import Read_BOSS
    >> data = Read_BOSS([arguments]).read()
    
where the instance of the Read_BOSS class is created using the [arguments] provided (see "Usage" section above for help on what these arguments are and what values they take by default).

### To run Read_BOSS from shell

To run Read_BOSS from shell simply type:

    $ python Read_BOSS.py [optional arguments]
    
To see the allowed optional arguments to be passed from the command line, type (or see "Usage" section):

    $ python Read_BOSS.py -h
 
 which outputs:
 
     usage: Read_BOSS.py [-h] [-gn GNAMENORTH] [-gs GNAMESOUTH] [-rn RNAMENORTH]
                         [-rs RNAMESOUTH] [-z ZLIMITS] [-l LENS] [-r SAVE_RDZWN]
                         [-x SAVE_XYZWN] [-s USE_SYS] [-hub HUBBLE]
                         [-m OMEGAMATTER]

     Code which reads the BOSS galaxy and random catalogs in fits format and
     returns text files containing RA, DEC, Z, W, N or X, Y, Z, W, N columns
     (allowing choice of cosmology). Reading is performed by the read() method
     defined in the Read_BOSS class.

     optional arguments:
       -h, --help            show this help message and exit
       -gn GNAMENORTH, --galnorth GNAMENORTH
                             name of Northern galaxy catalogue (with extension)
       -gs GNAMESOUTH, --galsouth GNAMESOUTH
                             name of Southern galaxy catalogue (with extension)
       -rn RNAMENORTH, --rannorth RNAMENORTH
                             name of Northern random catalogue (with extension)
       -rs RNAMESOUTH, --ransouth RNAMESOUTH
                             name of Southern random catalogue (with extension)
       -z ZLIMITS, --redshift ZLIMITS
                             upper and lower limit for redshift cut (given as a
                             list, default is [0.43, 0.7])
       -l LENS, --length LENS
                             number of galaxies and randoms after cut (given as a
                             list, default is [777202, 40096792])
       -r SAVE_RDZWN, --rdzwn SAVE_RDZWN
                             whether or not to save RA, DEC, Z, W, N files (answer
                             Yes, True', T, Y, 1 or No, False, F, N, 1)
       -x SAVE_XYZWN, --xyzwn SAVE_XYZWN
                             whether or not to save X, Y, Z, W, N files (answer
                             Yes, True', T, Y, 1 or No, False, F, N, 1)
       -s USE_SYS, --sys USE_SYS
                             whether or not to use systematics for the galaxy
                             catalogue (answer Yes, True', T, Y, 1 or No, False, F,
                             N, 1)
       -hub HUBBLE, --hubble HUBBLE
                             Hubble parameter, H0, in units of km/s/Mpc
       -m OMEGAMATTER, --matter OMEGAMATTER
                             Omega_matter

     For more information or help email sunny.vagnozzi@fysik.su.se.
     
### Examples
