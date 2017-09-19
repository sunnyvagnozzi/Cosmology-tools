# chi2_neutrino_hierarchy

A Python code which allows the combination of cosmological and neutrino oscillations data in order to determine the degree of preference for the normal mass hierarchy (NH) over the inverted mass hierarchy (IH). Specifically, the code combines posteriors on the sum of the three active neutrino masses obtained from cosmological data with the Chi^2 profiles for the solar and atmospheric mass splittings from NuFIT 3.0 (http://www.nu-fit.org/?q=node/139), and outputs the significance (provided as # of sigmas) at which NH is preferred over IH. The code provides the chi2_neutrino_hierarchy class, wherein the num_sigma() method calculates the desired preference for NH over IH.

## Usage

### To run the code from an ipython terminal or notebook

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

## Examples

Some examples on how to run the code from shell. If my M_nu posterior from cosmological data is in the Posteriors/ folder, and is called "planckTT_lowP_BAO_p_mnu.dat", plus I want the code to save the Chi^2 profiles (which will be saved in profile_chi2_nh.dat and profile_chi2_ih.dat), and finally I want the code to be chatty (i.e. to execute print statemements to keep me informed on what's going one), then I would run the code as follows:

     python chi2_neutrino_hierarchy.py -r Posteriors/planckTT_lowP_BAO -s Yes -v 1
     
If I want to use only oscillations data, and still want to save the Chi^2 profiles, but don't want the code to be chatty, I run:

     python chi2_neutrino_hierarchy.py -d 1

If you want to try out the code, in this same folder I have provided two cosmological posteriors you can use, "gsnpde_p_mnu.dat" and "lcdm_p_mnu.dat".

## Known issues and bugs

If the code is run with the flag "-d 2" (i.e. only cosmology data is used), the code will return a nan sigma preference for the NH. The reason for this is that cosmological data is not sensitive to the neutrino mass hierarchy (at least at present time and to leading order), and hence the cosmological posterior one uses is the same for both NH and IH. Therefore, the resulting \Delta \Chi ^2 = 0, which ends up being translated into a nan sigma preference for one hierarchy over the other (more properly, it should be 0 sigma). In a future release this will be fixed to account for a proper treatment of the preference for NH over IH with cosmology data alone, as per Hannestad & Schwetz, JCAP 11 (2016) 035, http://arxiv.org/abs/arXiv:1606.04691

## Things to keep in mind

When you run the code, make sure that you have included the following codes or files (which I provide in this folder) in the same folder in which you are running the code:

-*chisq2sigma.py*: a Python class used to convert a value of Delta Chi^2 to a significance in # of sigma, for any given number of degrees of freedom of the Chi^2 distribution

-*deltam2_21_nh.dat*: Chi^2 profiles for the solar mass splitting assuming NH provided by NuFIT 3.0

-*deltam2_32_nh.dat*: Chi^2 profiles for the atmospheric mass splitting assuming NH provided by NuFIT 3.0

-*deltam2_21_ih.dat*: Chi^2 profiles for the solar mass splitting assuming IH provided by NuFIT 3.0

-*deltam2_32_ih.dat*: Chi^2 profiles for the atmospheric mass splitting assuming IH provided by NuFIT 3.0

The code relies on standard Python libraries such as numpy, scipy, and argparse in order to run, so make sure you have those installed.
