# Driver program for creating a database of snapshots of the quasi-1D steady
# C-D nozzle flow problem, as specified by Nair & Balajewicz (2019)

import os
import numpy as np
import matplotlib.pylab as plt
from optparse import OptionParser

import Quasi1D_Steady as qs


parser = OptionParser(usage="usage: %prog -f filename -p pathOut -d debug")
parser.add_option('-p',dest='path',default='.',help='Output path')
parser.add_option('-v',dest='validate',action="store_true",default=False, \
    help='Validation (if supplied); else training')
(options, args) = parser.parse_args()


# Training database parameter set
mus_train = [0.5,0.875,1.25,1.625]   #Throat area values

# Parameter set and output file name based on whether 'Training' or 'Validation'
# database is to be generated
if not options.validate:    #Training
    outfile = os.path.join(options.path,'training_snapshots.npz')
    mus = mus_train
else:                       #Validation
    outfile = os.path.join(options.path,'testing_snapshots.npz')
    # Validation database is formed by placing two uniformly-distributed mu's
    # in every interval of mu's in the training database
    mus = []    # Validation database
    for im in range(len(mus_train)-1):
        mu_l = mus_train[im];     mu_r = mus_train[im+1]
        dmu = mu_r - mu_l
        mus.append(mu_l+dmu/3.)
        mus.append(mu_r-dmu/3.)

""" Quasi-1D steady shocked flow in C-D nozzle from Nair & Balajewicz 2019 """
xi = 0.         #Initial x-coordinates
xf = 10.        #Final x-coordinates
Am = 3          #C-D nozzle area at start (x = 0) and end (x = L) of x-domain
rhoin = 1.0     #Static density at inlet [kg/m^3]
pin_bar = 1.0   #Static pressure at inlet [bar]
pout_bar = 0.7  #Static pressure at outlet [bar]

bar2Pa = 1.     #Bar to Pascal conversion
pin = pin_bar*bar2Pa    #Static pressure at inlet [Pa]
pout = pout_bar*bar2Pa  #Static pressure at outlet [Pa]

Nx = 1001       # No. of x-grid points
xarr = np.linspace(xi,xf,Nx)    #Compose the x-grid

""" Allocate data structure for storing database """
nV = 3                              #No. of variables in solution
u_db = np.zeros((Nx,nV,len(mus)))   #Solution (all components) for all mu's
A_db = np.zeros((Nx,len(mus)))      #Area distribution for all mu's
# Obtain the problem-setup data structure (array) for the first parameter value
# just for knowing the size of the array
prblm_setup0 = qs.quasi_1D_steady_prblm_setup(rhoin,pin,pout,Am,mus[0], \
        NozzleKind=1)
# Allocate 2D array for storing problem setup data structures for all mu's
prblm_setups = np.zeros((len(prblm_setup0),len(mus)))

""" Calculate snapshots """
for im, mu in enumerate(mus):
    # Setup the problem
    prblm_setups[:,im] = qs.quasi_1D_steady_prblm_setup(rhoin,pin,pout,Am,mu, \
        NozzleKind=1)
    # Solve the quasi 1D steady flow problem in the above-defined C-D nozzle
    u_db[:,:,im] = qs.quasi_1D_steady_soln_in(prblm_setups[:,im],xarr)
    # Area distribution
    A_db[:,im] = qs.quasi_1D_steady_A_ovr_At(prblm_setups[:,im],xarr)*mu

""" Plot the snapshots """
legends = ['mu = '+str(mu) for mu in mus]   #Form the legends
qs.quasi_1D_steady_soln_plot(xarr,u_db,add_in=A_db,add_title='Area', \
    legends=legends)
plt.show()

# Save the solution
np.savez(outfile,prblm_setups=prblm_setups,xarr=xarr,A_db=A_db,u_db=u_db,mus=mus)
