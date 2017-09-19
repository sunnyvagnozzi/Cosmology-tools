# chi2_neutrino_hierarchy

A Python code which allows the combination of cosmological and neutrino oscillations data in order to determine the degree of preference for the normal mass hierarchy (NH) over the inverted mass hierarchy (IH). Specifically, the code combines posteriors on the sum of the three active neutrino masses obtained from cosmological data with the Chi^2 profiles for the solar and atmospheric mass splittings from NuFIT 3.0 (http://www.nu-fit.org/?q=node/139), and outputs the significance (provided as # of sigmas) at which NH is preferred over IH. The code provides the chi2_neutrino_hierarchy class, wherein the num_sigma() method calculates the desired preference for NH over IH.

## Usage

### Run the code from an ipython terminal or notebook

To run the code from an ipython terminal or notebook (in the same folder as chi2_neutrino_hierarchy.py), simply type:

    >> from chi2_neutrino_hierarchy import chi2_neutrino_hierarchy
    >> sigma = chi2_neutrino_hierarchy([optional arguments]).num_sigma()
    
where the instance of the chi2_neutrino_hierarchy class is created using the [optional arguments] provided (see the section below for help on what these arguments are and what values they take by default, or type "python chi2_neutrino_hierarchy.py -h" in the command line for more help).

### To run the code from shell

To run the code from shell simply type:

    $ python chi2_neutrino_hierarchy.py [optional arguments]
    
To see the allowed optional arguments to be passed from the command line, type:

    $ python chi2_neutrino_hierarchy.py -h
 
 which outputs:
 
     usage: chi2_neutrino_hierarchy.py [-h] [-r ROOT] [-s SAVE_PROFILES]
                                       [-d DATASET] [-v VERBOSE]

     ** Code which combines cosmology and neutrino oscillation data from NuFIT 3.0
     to determine the # of sigmas at which the normal hierarchy is preferred over
     the inverted hierarchy. **

     optional arguments:
       -h, --help            show this help message and exit
       -r ROOT, --root ROOT  Root of the M_nu posterior file ([root]_p_mnu.dat),
                             includes the path here if required (default:
                             'gsnpde').
       -s SAVE_PROFILES, --save SAVE_PROFILES
                             Whether or not to save Chi^2 profiles against M_nu
                             (answer Yes, True, T, Y, 1 or No, False, F, N, 1,
                             default is False)
       -d DATASET, --dataset DATASET
                             Which datasets to use to calculate the Chi^2 0:
                             oscillations and cosmology (default) 1: only
                             oscillations 2: only cosmology
       -v VERBOSE, --verbose VERBOSE
                             Controls chattiness of code: 0 (default) for non
                             chatty, >0 for chatty.

     ** Developed by Sunny Vagnozzi, 2017. For more information, help, or for
     reporting bugs please contact sunny.vagnozzi@fysik.su.se. For help with the
     usage of the code, type in the command line: python chi2_neutrino_hierarchy.py
     -h Known bugs: when running with the flag "-d 2" (i.e. when only cosmological
     data is usecd) the code will return nan sigma since cosmological data is not
     sensitive to the neutrino mass ordering and hence Delta Chi^2 (NH-IH) = 0.
     Will be improved in a future release through the treatment of Hannestad &
     Schwetz 2016. **
