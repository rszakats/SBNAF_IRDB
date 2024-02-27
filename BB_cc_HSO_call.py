# Python version of the Black Body colour correction code.
#
# - Victor Ali-Lagoa, MPE, 2018-01
# RS - mod to be callable finction from python3
import sys
'''
###################################################
if len(sys.argv) < 2:
    print "
Calculates the black-body colour correction for a given temperature
Run as:
> python {0} temperature_in_K
Output:
T[K]  K(12.0) K(25.0) K(60.0) K(100.)""".format(sys.argv[0])"
    sys.exit(1)
'''
import math
import numpy as np

###################################################
# CONSTANTS, PARAMETERS, etc.

# Reference wavelengths (not necessarily isophotal wavelengths)
lambda_ref=[70.0, 100.0, 160.0]

# Phys. constants for Black Body SED
h_planck = float(6.62606896e-34)  # in J.s
c_light  = float(2.997925e+8)     # in m/s
K_B      = float(1.3806504e-23)   # in J/K

cB = 2.0*1.0e-10*h_planck*c_light*c_light
# with this constant,
# the B(lambda, T) will be in W/cm^2/micron, i.e.,
# the units they use in Wright et al. (2010)
# because:
#
# J.s x m/s x m/s = J/s m m = W.100 cm.100 cm/lambda5 =
# = W 100 cm 100 cm / (100 cm * 100 cm * 100 cm * 100 cm * 1000000 micron) =
# = W/cm2/1e10micron = 1e-10 W/cm2/micron

###################################################
# VARIABLES
# one index per filter:
colour_correction     = np.ones(4)
num_rows = np.zeros(4, dtype = np.int64)
l_filt  = np.zeros((4,3000)) - 1.0
r_filt  = np.zeros((4,3000))
re_filt = np.zeros((4,3000)) - 1.0
f_data  = np.zeros((4,3000)) - 1.0

###################################################
# FUNCTIONS
###################################################
def interp(x,x1,x2,y1,y2):
    slope=(float(y2)-float(y1))/(float(x2)-float(x1))
    return float(y1) + (float(x)-float(x1))*slope

filenames=['PACS_B.dat','PACS_G.dat','PACS_R.dat']
def read_filters():
    # read the
    # wavelengths and responses for all the filters
    for K in range(3):
        filt_file = filenames[K]
        filt_data = np.genfromtxt(filt_file, autostrip = True, comments = "#")
        n_rows = len(filt_data)

        l_read  = filt_data[:,0]
        r_read  = filt_data[:,1]
        # re_read = filt_data[:,2]
        # no error column in the interpolated version of the filters
        num_rows[K] = n_rows

        for J in range(n_rows):
            l_filt[K,J]  = l_read[J]
            r_filt[K,J]  = r_read[J]
            # re_filt[K,J] = re_read[J]
            # no error column in the interpolated version of the filters

def colour_correct(K,n_rows,lambda_,response_,filt_wvl,T_in):

    c_in = float((1.0e+6)*h_planck*c_light/(K_B*T_in))   #in W/cm^2/micron
    # copy to define the black body spectrum and reference spectrum arrays:
    B_lambda_ = lambda_.copy()
    ref_      = lambda_.copy()

    # Prepare arrays:
    lambda_5 = lambda_**5
    inv_lambda_5 = (1.0e+30)/lambda_5
    # in m^(-5)

    # Fill in the B_lambda_ and ref_ arrays
    for m in range(len(lambda_)):
        frac1 = 1.0/lambda_[m]
        frac5 = inv_lambda_5[m]

        B_lambda_[m] = frac5*cB*(1.0/(math.exp(c_in*frac1) - 1.0))
        ref_[m] = frac1

    # We find the indices corresponding to the wavelengths enclosing
    # the reference wavelength of the filter filt_wvl
    for ii in range(len(lambda_)-1):
        if lambda_[ii] <= filt_wvl and lambda_[ii+1] >= filt_wvl: break

    # and then we interpolate the input and reference SEDs at filt_wvl
    B_filt_wvl = interp(filt_wvl,
                        lambda_[ii],lambda_[ii+1],
                        B_lambda_[ii],B_lambda_[ii+1])
    ref_filt_wvl= interp(filt_wvl,
                         lambda_[ii],lambda_[ii+1],
                         ref_[ii],ref_[ii+1])

    # Now we normalize B_lambda_ and ref_ to be 1 at filt_wvl:
    B_lambda_norm = B_lambda_.copy()
    ref_norm = B_lambda_.copy()
    for i in range(len(B_lambda_)):
        B_lambda_norm[i]  = (B_lambda_[i]/B_filt_wvl)
        ref_norm[i] = (ref_[i]/ref_filt_wvl)

    #Finally, we integrate to find in-band fluxes
    band_B_lambda_ = np.trapz(response_*B_lambda_norm,lambda_)
    band_ref_ = np.trapz(response_*ref_norm,lambda_)

    #print("{0:6.4f} ".format(band_B_lambda_/band_ref_))
    return "{0:6.4f} ".format(band_B_lambda_/band_ref_)

def ccorrect_hso(T_in):

    # read the wavelengths and responses  (l_filt and r_filt)
    # and store the number of rows (num_rows) for each filter
    read_filters()
    cc=[]
    #print("{0:<4.0f}".format(T_in))
    for K in range(0,3):
        cc.append(colour_correct(K,num_rows[K],l_filt[K,0:num_rows[K]],r_filt[K,0:num_rows[K]],lambda_ref[K],T_in))
    return(cc)
