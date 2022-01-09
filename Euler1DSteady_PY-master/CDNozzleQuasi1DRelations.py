# Quasi-1D flow relations for converging-diverging nozzles

import numpy as np

from MyPythonCodes.tools.fzeros import fzero_bisect
from MyPythonCodes.CompressibleFlow import IsentropicFlow as ifr
from MyPythonCodes.CompressibleFlow import NormalShock as nsr

gam_dflt = 1.4          #Default ratio of specific heats

# Finds ratio of local area to throat area (A/At) from local Mach no. using eqn.
# 5.20 of Anderson (2004)
def AratioFromM(M,gam=gam_dflt):
    Aratio = np.sqrt(pow(2./(gam+1)*(1+(gam-1)/2.*M**2),(gam+1)/(gam-1))/M**2)
    return Aratio

# Inverts the above function; i.e., finds Mach no. from area ratio; Fig 5.13 of
# Anderson (2004) shows that there are 2 solutions; the branch that is picked
# depends on whether the flag is 'subsonic' (default) or 'supersonic'.
# For safety, we use the bisection search method.
def MFromAratio(Aratio,flag='subsonic',gam=gam_dflt):
    if flag.lower() == 'subsonic':
        M_l = 1.e-6 #Left bracketing solution
        M_r = 1.    #Right bracketing solution
    elif flag.lower() == 'supersonic':
        M_r = 25    #Left bracketing solution
        M_l = 1.    #Right bracketing solution
    else:
        raise Exception('flag is neither subsonic nor supersonic')
    # Define function whose root is to be found
    func = lambda M : AratioFromM(M,gam=gam) - Aratio
    # Find root of above function by bisection method
    M, info, ier, mesg = fzero_bisect(func,M_l,M_r,full_output=1)
    if ier != 1:
        raise Exception(mesg)
    return M

# Given the ratio of exit area to throat area (Ae), returns back pressure (as a
# fraction of upstream stagnation pressure) needed for sonic flow at throat and
# SUBSONIC flow at exit. It uses MFromAratio to get exit Mach no. from exit area
# ratio; then it uses isentropic relation for static pressure as a fraction of
# stagnation pressure for given Mach no. - i.e. eqn. 3.30 of Anderson (2004)
def calc_pb_subsonic(Ae,gam=gam_dflt):
    return 1./ifr.IF_p(MFromAratio(Ae,'subsonic',gam=gam))
    
# Same as above, but for SUPERSONIC flow at exit
def calc_pb_supersonic(Ae,gam=gam_dflt):
    return 1./ifr.IF_p(MFromAratio(Ae,'supersonic',gam=gam))

# Given the area ratio at the exit, returns back pressure needed for sonic
# flow at throat and normal shock standing at exit
def calc_pb_normalshockexit(Ae,gam=gam_dflt):
    pe1 = calc_pb_supersonic(Ae,gam=gam) #Static pressure upstream of shock
    Me1 = MFromAratio(Ae,'supersonic',gam=gam) #Mach no. upstream of shock
    return pe1*nsr.NS_p2p1(Me1,gam=gam) #Static pressure downstream of shock

# Does calculations for a shock standing in the diverging section of a nozzle.
# Given the nozzle area at the shock and the exit area, both referred to the
# pre-shock sonic flow area (i.e., Astar1), returns the post-shock sonic flow
# area and various flow properties at exit
def calc_for_shock_in_divSec(As_ovr_Astar1,Ae_ovr_Astar1,gam=gam_dflt):
    Ms1 = MFromAratio(As_ovr_Astar1,'supersonic',gam=gam) #Mach no. before shock
    Ms2 = nsr.NS_M2(Ms1,gam=gam) #Mach no. downstream of shock
    As_ovr_Astar2 = AratioFromM(Ms2,gam=gam) #Area ratio downstream of shock
    Astar2_ovr_Astar1 = As_ovr_Astar1/As_ovr_Astar2 #Sonic area ratios across shock
    Ae_ovr_Astar2 = Ae_ovr_Astar1/Astar2_ovr_Astar1 #Exit area ratio
    p02 = nsr.NS_po2po1(Ms1,gam=gam) #Total pressure downstream of shock
    Me = MFromAratio(Ae_ovr_Astar2,'subsonic',gam=gam) #Mach no. at exit
    pe = p02/ifr.IF_p(Me,gam=gam) #Static pressure at exit
    return Astar2_ovr_Astar1, p02, Me, pe
