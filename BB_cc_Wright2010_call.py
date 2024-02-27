# Python version of the Black Body colour correction code.
#
# - Victor Ali-Lagoa, MPE, 2016-12

# import sys
# import subprocess
import math
import numpy as np
# from StringIO import StringIO
# import scipy as sp
# import matplotlib.pyplot as plt


# Red-blue calibrator discrepancy isoph. wavelengths
# cheat_iso_lambda=[3.3526,  4.6028, 11.0984, 22.6405]

# Ground isoph. wavelengths (Wright et al. 2010)
cheat_iso_lambda = [3.3526, 4.6028, 11.5608, 22.0883]


# CONSTANTS AND PARAMETERS
###################################################
# Phys. constants
h_planck = float(6.62606896e-34)  # in J.s
c_light  = float(2.997925e+8)     # in m/s
K_B      = float(1.3806504e-23)   # in J/K

# Vega spectrum constants (see reference in Wright et al. )
c1 = float(1.0158000e-16)
c2 = float(0.0083000)
c3 = float(8.8910000)             # in microns
# temperature for the Vega continuum
# (2) Wright et al. (2010)
T_V = float(14454.000)                             #in Kelvin
c_V = float((1.0e+6)*h_planck*c_light/(K_B*T_V))   #in W/cm^2/micron

cB = 2.0*1.0e-10*h_planck*c_light*c_light
# with this constant,
# the B(lambda, T) will be in W/cm^2/micron, i.e.,
# the units they use in Wright et al. (2010)
# because:
#
# J.s x m/s x m/s = J/s m m = W.100 cm.100 cm/lambda5 =
# = W 100 cm 100 cm / (100 cm * 100 cm * 100 cm * 100 cm * 1000000 micron) =
# = W/cm2/1e10micron = 1e-10 W/cm2/micron

# VARIABLES
###################################################
# one index for each filter:
passband_flux_V   = np.zeros(4)
passband_constant = np.zeros(4)
flux_V            = np.zeros(4)
flux_iso_V        = np.zeros(4)
iso_index         = np.zeros(4)
lambda_iso_V      = np.zeros(4)
passband_flux     = np.zeros(4)
mean_flux         = np.zeros(4)
flux_iso_tm       = np.zeros(4)

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

def read_filters():
    # read the
    # wavelengths and responses for all the filters
    filenames=["R_lambda_W1.txt", "R_lambda_W2.txt", "R_lambda_W3.txt","R_lambda_W4.txt"]
    for K in range(len(filenames)):
        # print("readfilters K:", K)
        filt_file = filenames[K]
        filt_data = np.genfromtxt(filt_file, autostrip = True, comments = "#")
        n_rows = len(filt_data)

        l_read  = filt_data[:,0]
        r_read  = filt_data[:,1]
#        re_read = filt_data[:,2]
# no error column in the interpolated version of the filters
        num_rows[K] = n_rows

        for J in range(n_rows):
            l_filt[K,J]  = l_read[J]
            r_filt[K,J]  = r_read[J]
#            re_filt[K,J] = re_read[J]
# no error column in the interpolated version of the filters

def Vega_isoph_fluxes(K,n_rows,lambda_,response_,l_iso_V):

    lambda_5 = lambda_**5
    inv_lambda_5 = (1.0e+30)/lambda_5
    # in m^(-5)

    B_lambda_V = lambda_.copy()
    flux_V     = lambda_.copy()
    for m in (len(lambda_)):
        frac5 = inv_lambda_5[m]
        frac1 = 1.0/lambda_[m]
        B_lambda_V[m] = frac5*cB*(1.0/(math.exp(c_V*frac1) - 1.0))
        flux_V[m] = c1*B_lambda_V[m]*(
            1.0 - c2*math.log(lambda_[m]/c3)*math.log(lambda_[m]/c3)
        )
#        print "xV{0:1d} {1:6.4f} {2:10.7e}".format(
#            K+1,lambda_[m], flux_V[m])

    # We check the isophotal fluxes here:
    # passband flux of Vega in the K filter
    # and denom = passband_constant[K]
    num   = sum(response_*lambda_*flux_V)
    denom = sum(response_*lambda_)
    in_band_V = num/denom

    # To use the input lambda_iso,
    # (for the old version, see commented lines at the end)
    l_iso_V=cheat_iso_lambda[K]

    # we need the indices that correspond to the wavelengths enclosing
    # the value of lambda_iso:
    for ii in range(len(lambda_)-1):
        if lambda_[ii] <= l_iso_V and lambda_[ii+1] >= l_iso_V: break

    # so we interpolate
    f_iso_V=interp(l_iso_V,
                   lambda_[ii],lambda_[ii+1],
                   flux_V[ii],flux_V[ii+1])

    # To compare with Wright et al., we need the conversion factor
    convfact=l_iso_V*l_iso_V/c_light*1e24
    print_f_iso_V=f_iso_V*convfact
    if (K==3):
        # for band W4 (K==3), we add a 2.7 per cent
        # (see Wright et al. 2010)
        print_f_iso_V*=1.027

    '''
    print "W{0}: {1: 8.4f} F_iso(W{0}) = {2:8.4f} ".format(
        K+1, l_iso_V, print_f_iso_V),
    print "lii={0:5.2f} lii_1={1:5.2f} ".format(lambda_[ii],
                                                lambda_[ii+1]),
    print "fii={0:9.4f} fii_1={1:9.4f} inband={2:9.4f}".format(
        flux_V[ii]*convfact, flux_V[ii+1]*convfact, in_band_V*convfact)

    print "W{0:1d} SUM(lambda VegaNorm_lambda)D_lambda = {1:10.7e}  Lambda_iso = {2:6.4f}".format(K+1, num/f_iso_V,l_iso_V)
    '''

def colour_correct_2(K,n_rows,lambda_,response_,l_iso_V,flag_1,flag_2):

    lambda_5 = lambda_**5
    inv_lambda_5 = (1.0e+30)/lambda_5
    # in m^(-5)

    # Black body
    B_lambda_ = lambda_.copy()
    flux_     = lambda_.copy()

    # Vega
    B_lambda_V = lambda_.copy()
    flux_V     = lambda_.copy()
    for m in range(len(lambda_)):
        frac5 = inv_lambda_5[m]
        frac1 = 1.0/lambda_[m]

        # Vega
        B_lambda_V[m] = frac5*cB*(1.0/(math.exp(c_V*frac1) - 1.0))
        flux_V[m] = c1*B_lambda_V[m]*(
            1.0 - c2*math.log(lambda_[m]/c3)*math.log(lambda_[m]/c3)
        )

        # Input SED (black body or other)
        B_lambda_[m] = frac5*cB*(1.0/(math.exp(c_in*frac1) - 1.0))
        ##  lambda**2 = nu**(-4)
        # B_lambda_[m] = lambda_[m]**2
        ## constant lambda = nu*(-2)
        # B_lambda_[m] = 1
        ## lambda**(-2) = constant nu
        # B_lambda_[m] = 1./(lambda_[m]*lambda_[m])
        ## lambda**(-4) = nu**2
        # B_lambda_[m] = 1./(lambda_[m]*lambda_[m]*lambda_[m]*lambda_[m])

    l_iso_V=cheat_iso_lambda[K]
    # We find the indices corresponding to the wavelengths enclosing l_iso_V
    for ii in range(len(lambda_)-1):
        if lambda_[ii] <= l_iso_V and lambda_[ii+1] >= l_iso_V: break

    # We interpolate the Vega spectrum at l_iso_V
    f_iso_V=interp(l_iso_V,
                   lambda_[ii],lambda_[ii+1],
                   flux_V[ii],flux_V[ii+1])

    # passband flux of Vega in the K filter
    # normalised to 1 at l_iso_V
    num_V   = sum(response_*lambda_*flux_V)
    norm_in_band_V = num_V/f_iso_V

    # We interpolate input SED at l_iso_V:
    tmf_li= interp(l_iso_V,
                   lambda_[ii],lambda_[ii+1],
                   B_lambda_[ii],B_lambda_[ii+1])

    # normalize B_lambda to be 1 at lambda_iso:
    tm_flux_norm  = B_lambda_.copy()
    for i in range(len(B_lambda_)):
        tm_flux_norm[i]  = (B_lambda_[i]/tmf_li)

    # and integrate to find in-band flux of input SED normalised at l_iso_V
    norm_in_band_flux_l = sum(response_*lambda_*tm_flux_norm)

    '''
    if flag_1:
        print "W{0:1d} SUM(lambda    FNorm_lambda)D_lambda = {1:10.7e}".format(
            K+1, norm_in_band_flux_l)
    '''
    colour_correction[K] = norm_in_band_flux_l/norm_in_band_V



def colour_correct(K,n_rows,lambda_,response_,l_iso_V,T_in):
    c_in = float((1.0e+6)*h_planck*c_light/(K_B*T_in))   #in W/cm^2/micron
    lambda_5 = lambda_**5
    inv_lambda_5 = (1.0e+30)/lambda_5
    # in m^(-5)

    B_lambda_ = lambda_.copy()
    flux_     = lambda_.copy()
    for m in range(len(lambda_)):
        frac5 = inv_lambda_5[m]
        frac1 = 1.0/lambda_[m]
        B_lambda_[m] = frac5*cB*(1.0/(math.exp(c_in*frac1) - 1.0))
        ##  lambda**2 = nu**(-4)
        # B_lambda_[m] = lambda_[m]**2
        ## constant lambda = nu*(-2)
        # B_lambda_[m] = 1
        ## lambda**(-2) = constant nu
        # B_lambda_[m] = 1./(lambda_[m]*lambda_[m])
        ## lambda**(-4) = nu**2
        # B_lambda_[m] = 1./(lambda_[m]*lambda_[m]*lambda_[m]*lambda_[m])

#        print "xF{0:1d} {1:6.4f} {2:10.7e} {3:10.7e}".format(
#                K, lambda_[m], B_lambda_[m],
#                B_lambda_[m]*lambda_[m]*response_[m])
#        if flag_2:
#            print "xxF{0:1d} {1:6.4f} {2:10.7e} {3:10.7e}".format(
#                K, lambda_[m], B_lambda_[m],
#                B_lambda_[m]*lambda_[m]*response_[m])

    l_iso_V=cheat_iso_lambda[K]
    # We find the indices corresponding to the wavelengths enclosing l_iso_V
    for ii in range(len(lambda_)-1):
        if lambda_[ii] <= l_iso_V and lambda_[ii+1] >= l_iso_V: break
    # We interpolate input SED at l_iso_V:
    tmf_li= interp(l_iso_V,
                   lambda_[ii],lambda_[ii+1],
                   B_lambda_[ii],B_lambda_[ii+1])

    # normalize B_lambda to be 1 at lambda_iso:
    tm_flux_norm  = B_lambda_.copy()
    for i in range(len(B_lambda_)):
        tm_flux_norm[i]  = (B_lambda_[i]/tmf_li)

    # and integrate to find in-band flux of input SED normalised at l_iso_V
    band_tm_flux_l = sum(response_*lambda_*tm_flux_norm)

    # The band's Delta_Lambda*R_Lambda*d_Lambda is the denominator:
    denom = sum(response_*lambda_)
    '''
    if flag_1:
        print "W{0:1d} SUM(lambda    FNorm_lambda)D_lambda = {1:10.7e}".format(
            K+1, band_tm_flux_l)
    '''
    #colour_correction[K] = band_tm_flux_l/(denom)
    return "{0:6.6f} ".format(band_tm_flux_l/(denom))


flaG1  = True
flagG2 = True

#T_in = float(sys.argv[1])


# read de wavelengths and responses  (l_filt and r_filt)
# and store the number of rows (num_rows) for each filter


def ccorrect_wise(T_in):
    read_filters()
    cc = []
    for K in range(len(cheat_iso_lambda)):
        # print("K:",K)
        '''
        cc.append(Vega_isoph_fluxes(K,
                          num_rows[K],
                          l_filt[K,0:num_rows[K]],
                          r_filt[K,0:num_rows[K]],
                          cheat_iso_lambda[K]))
        '''
        cc.append(colour_correct(K,
                      num_rows[K],
                      l_filt[K,0:num_rows[K]],
                      r_filt[K,0:num_rows[K]],
                      cheat_iso_lambda[K],T_in))
    return(cc)



## Old estimate of lambda_iso based on
## Definition of isophotal flux: where flux_V(lambda)=flux_iso_V
# flux_iso_V[K] = loc_fi
#for ii in range(len(lambda_)-1):
#    if loc_fi <= flux_V[ii] and loc_fi >= flux_V[ii+1]: break
## got the index where the above contition is met: ii
##print "old ii", ii, lambda_[ii]
