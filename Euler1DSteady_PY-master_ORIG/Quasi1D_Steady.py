#Solver for quasi 1-D steady flow in a C-D nozzle with possible shock

import numpy as np
import matplotlib.pylab as plt

from MyPythonCodes.tools import fzero_bisect
from MyPythonCodes.CompressibleFlow import IsentropicFlow as ifr
from MyPythonCodes.tools import Findiff_Taylor_uniform

import CDNozzleQuasi1DRelations as cdn
import Quasi1D_Steady_Nozzles

gam_dflt = 1.4          #Default ratio of specific heats

""" 
Dict mapping 'NozzleKind' integer to corresponding 'Quasi1D_Steady_Nozzles' 
function names 
"""
quasi_1D_steady_funs_A_ovr_At = {
    1 : 'fun_A_over_At_parabolic',
    }

def quasi_1D_steady_prblm_setup(rhoin,pin,pout,Am,At,NozzleKind=1,gam=gam_dflt):
    """
    Setup the problem - compose an array that encodes the various parameters

    INPUTS:
    rhoin       : Inlet density
    pin         : Inlet pressure
    pout        : Outlet pressure
    Am          : Inlet area
    At          : Throat area
    NozzleKind  : Kind of nozzle
    gam         : Sp. heat ratio
    
    OUTPUTS:
    prblm_setup : Array consisting of above parameters in the same order
    """
    return np.array([rhoin,pin,pout,Am,At,NozzleKind,gam])
#enddef quasi_1D_steady_prblm_setup


def quasi_1D_steady_fun_A_ovr_At(prblm_setup,xi,xf):
    """
    Obtain handle of function that will yield the local to throat area ratio of
    the C-D nozzle for input grid coordinates

    INPUTS:
    prblm_setup : Output of 'quasi_1D_steady_prblm_setup' function call
    xi          : Inlet grid coordinate
    xf          : Outlet grid coordinate
    
    OUTPUTS:
    fun_A_ovr_At : Handle to function giving area ratio for input grid 
                   coordinates
    """
    # Parse the problem setup
    Am          = prblm_setup[3]        #Inlet area
    At          = prblm_setup[4]        #Throat area
    NozzleKind  = int(prblm_setup[5])   #Kind of nozzle
    # Get function handle that will define the function for obtaining A/At
    def_fun_A_ovr_At = getattr(Quasi1D_Steady_Nozzles, \
        quasi_1D_steady_funs_A_ovr_At[NozzleKind])
    # Define the function for obtaining A/At; this needs the first & last grid
    # coordinates, and the inlet & throat areas
    return def_fun_A_ovr_At(xi,xf,Am,At)
#enddef quasi_1D_steady_fun_A_ovr_At


def quasi_1D_steady_A_ovr_At(prblm_setup,xarr):
    """
    Determine the local to throat area ratio of the C-D nozzle for input grid
    coordinates

    INPUTS:
    prblm_setup : Output of 'quasi_1D_steady_prblm_setup' function call
    xarr        : Array of grid coordinates along axis of nozzle
    
    OUTPUTS:
    Aarr : Area ratio array corresponding to 'xarr'
    """
    fun_A_ovr_At = quasi_1D_steady_fun_A_ovr_At(prblm_setup,xarr[0],xarr[-1])
    return fun_A_ovr_At(xarr)
#enddef quasi_1D_steady_A_ovr_At


def quasi_1D_steady_soln_basic(pb,xarr,fun_A_ovr_At,gam=gam_dflt):
    """
    Basic solver for quasi 1-D steady flow in a C-D nozzle w/ possible shock
    
    INPUTS:
    pb           : Ratio of exit (back) pressure to inlet stagnation pressure
    xarr         : Array of grid coordinates along axis of nozzle
    fun_A_ovr_At : Function giving ratio of local to throat area for any 'x'
    gam          : Ratio of specific heats (default is 'gam_dflt')
    
    OUTPUTS:
    soln : 2D array with x-grid points along 0th dimension and the following
           solution variables along 1st dimension
           i)   Local Mach no. for each entry of 'xarr'
           ii)  Ratio of local pressure to inlet stagnation pressure for 'xarr'
           iii) Ratio of local temperature to inlet stagnation temperature for
                'xarr'
    """
    Nx = len(xarr)
    soln = np.zeros((Nx,3))

    ####
    # In the following, we refer to A/A* as the area ratio exclusively.
    # If the throat isn't choked, then A/At != A/A*
    ####
    
    A_ovr_At = fun_A_ovr_At(xarr) #Ratio of cross-sectional area to throat area
    IAt = np.argmin(A_ovr_At) #Index of throat section
    Ae_ovr_At = A_ovr_At[-1]  #Ratio of exit area to throat area
    pb_subsonic = cdn.calc_pb_subsonic(Ae_ovr_At,gam=gam)
    pb_normalshockexit = cdn.calc_pb_normalshockexit(Ae_ovr_At,gam=gam)
    
    if pb >= pb_subsonic:
        """ Subsonic flow throughout (throat is not choked; At != A*) """
        Me = ifr.IF_unp(1./pb,gam=gam) #Exit Mach no. (necessarily subsonic)
        if (Me > 1.):
            raise Exception('For pb >= pb_subsonic (= '+str(pb_subsonic) \
                +'), Me (= '+str(Me)+') > 1!')
        Ae_ovr_Astar = cdn.AratioFromM(Me,gam=gam) #Exit area ratio (Ae/A*)
        if (Ae_ovr_Astar < Ae_ovr_At):
            raise Exception('For pb >= pb_subsonic (= '+str(pb_subsonic) \
                +'), At < Astar!')
        A_ovr_Astar = A_ovr_At*(Ae_ovr_Astar/Ae_ovr_At) #Area ratios
        for ix in range(Nx):  #Go thru each axial station to find local Mach no.
            # Actually find subsonic Mach no. for local area ratio (A/A*)
            soln[ix,0] = cdn.MFromAratio(A_ovr_Astar[ix],'subsonic',gam=gam)
        soln[:,1] = 1./ifr.IF_p(soln[:,0],gam=gam) #Corresponding pressure
    elif pb > pb_normalshockexit:
        """ Normal shock exists somewhere within the diverging section """
        ### We first find the location of the shock, xs, by using bisection
        # Initial 'left' and 'right' guesses are throat & exit, respectively
        # We know that pb is bracketed above and below by the static pressures
        # corresponding to the 'right' and 'left' guesses of 'xs', respectively
        xl = xarr[IAt];     xr = xarr[-1];
        # Define function whose root is to be found, such that if shock location
        # is x=xs, then exit pressure pe = pb, the specified back pressure; this
        # function returns the difference pe - pb for a given xs
        def func(xs):
            # Determine area ratio at current guess of shock location. N.B. this
            # refers to sonic area of pre-shock (and not post-shock) flow.
            As_ovr_Astar1 = fun_A_ovr_At(xs) # Area ratio (upstream of) at shock
            # Find pe/p01, p02/p01 and A2*/A1* for area ratio at current guess
            # shock location and exit area, both referred to actual throat area
            Astar2_ovr_Astar1, p02, Me, pe \
                = cdn.calc_for_shock_in_divSec(As_ovr_Astar1,Ae_ovr_At,gam=gam)
            return pe - pb
        # Find root of above function by bisection method
        xs, info, ier, mesg = fzero_bisect(func,xl,xr,full_output=1)
        if ier != 1:    #Unsuccessful exit from root-finding method
            raise Exception(mesg)
        # Determine area ratio at shock location. N.B. this refers to sonic area
        # of pre-shock (and not post-shock) flow.
        As_ovr_Astar1 = fun_A_ovr_At(xs) # Area ratio (upstream of) at shock
        # Find pe/p01, p02/p01 and A2*/A1* for area at shock location and exit
        # area, both referred to actual throat area
        Astar2_ovr_Astar1, p02, Me, pe \
            = cdn.calc_for_shock_in_divSec(As_ovr_Astar1,Ae_ovr_At,gam=gam)
        # Index of point just to the left of (or on) the shock
        isl = np.argmin(abs(xarr - xs))
        if xarr[isl] > xs:
            isl = isl - 1
        ### Calculate Mach no. & pressure in the converging and the unshocked
        ### diverging sections
        for ix in range(isl+1):
            if ix <= IAt:   #In converging section
                soln[ix,0] = cdn.MFromAratio(A_ovr_At[ix],'subsonic',gam=gam)
            else:   #In diverging section before shock
                soln[ix,0] = cdn.MFromAratio(A_ovr_At[ix],'supersonic',gam=gam)
        # Calculate ratio of corresponding static to stagnation pressures
        soln[:isl+1,1] = 1./ifr.IF_p(soln[:isl+1,0],gam=gam)
        # Calculate Mach no. in diverging section downstream of the shock
        # (N.B.: A* has changed due to the non-isentropic shock; account for it)
        for ix in range(isl+1,Nx):
            soln[ix,0] = cdn.MFromAratio(A_ovr_At[ix]/Astar2_ovr_Astar1, \
                'subsonic',gam=gam)
        # Calculate corresponding static pressure, using 'p02'
        soln[isl+1:,1] = p02/ifr.IF_p(soln[isl+1:,0],gam=gam)
    else:
        # Flow is isentropic throughout the duct, with possible non-isentropy
        # at exit plane or outside
        for ix in range(Nx):
            if ix <= IAt:   #In converging section
                soln[ix,0] = cdn.MFromAratio(A_ovr_At[ix],'subsonic',gam=gam)
            else:   #In diverging section before shock
                soln[ix,0] = cdn.MFromAratio(A_ovr_At[ix],'supersonic',gam=gam)
        soln[:,1] = 1./ifr.IF_p(soln[:,0],gam=gam) #Corresponding p/p0
    # Stagnation temperature is constant throughout the nozzle (it doesn't
    # change across a shock); so we can directly calculate ratio of static to
    # INLET stagnation temperatures everywhere from the known local Mach nos. 
    soln[:,2] = 1./ifr.IF_T(soln[:,0],gam=gam)
    return soln
#enddef quasi_1D_steady_soln_basic


def quasi_1D_steady_soln_in(prblm_setup,xarr):
    """
    Solver for quasi 1-D steady flow in a C-D nozzle w/ possible shock; builds
    off of 'quasi_1D_steady_soln_basic', with modified inputs (inlet static
    conditions) and outputs (conserved flow variables)
    
    INPUTS:
    prblm_setup : Output of 'quasi_1D_steady_prblm_setup' function call
    xarr        : Array of grid coordinates along axis of nozzle
    
    OUTPUTS:
    soln : 2D array with x-grid points along 0th dimension and the conserved
           flow variables along 1st dimension
           i)   rho (density) for each entry of 'xarr'
           ii)  rho*u (momentum) for each entry of 'xarr'
           iii) rho*E (energy) for each entry of 'xarr'
    """
    # Parse the problem setup
    rhoin       = prblm_setup[0]        #Inlet density
    pin         = prblm_setup[1]        #Inlet pressure
    pout        = prblm_setup[2]        #Outlet pressure
    gam         = prblm_setup[6]        #Sp. heat ratio
    # Obtain handle of function that will yield the local to throat area ratio
    # for input grid coordinates; this needs the problem-setup data structure
    # as well as the inlet and outlet grid coordinates
    # coordinates, and the inlet & throat areas
    fun_A_ovr_At = quasi_1D_steady_fun_A_ovr_At(prblm_setup,xarr[0],xarr[-1])
    # Obtain stagnation pressure at inlet (needed by quasi_1D_steady_soln_basic)
    Arr = fun_A_ovr_At(xarr)    #Area ratio for all x-grid points
    A_in = Arr[0]               #Inlet area
    A_t = np.min(Arr)           #Throat area
    Min = cdn.MFromAratio(A_in/A_t,flag='subsonic',gam=gam)    #Inlet Mach no.
    p0in = pin*ifr.IF_p(Min) #Inlet stagnation pressure
    pout_by_p0in = pout/p0in  #Outlet static to inlet stagnation pressure
    # Get solution in default variables: Mach no., static to inlet stagnation
    # pressure, and static to inlet stagnation temperature
    soln_0 = quasi_1D_steady_soln_basic(pout_by_p0in,xarr,fun_A_ovr_At,gam=gam)
    # Retrieve solution elements (1D arrays) from 2D array
    Ms = soln_0[:,0];   ps_by_p0in = soln_0[:,1];   Ts_by_T0in = soln_0[:,2]
    # Pre-allocate storage for output solution array
    soln = np.zeros_like(soln_0)
    # 0th column of solution is
    #   rho = (rho/rho_0_in)*(rho0in/rhoin)*rhoin
    #       = (p/p0in)/(T/T0_in)*(rho0in/rhoin)*rhoin
    soln[:,0] = ps_by_p0in/Ts_by_T0in*ifr.IF_d(Min,gam=gam)*rhoin
    # 1st column of solution is
    #   rho*u = rho*M*a = rho*M*sqrt(gam*p/rho) = M*sqrt(gam*p*rho)
    #         = M*sqrt(gam*p0in*rho*(p/p0in))
    soln[:,1] = Ms*np.sqrt(gam*p0in*soln[:,0]*ps_by_p0in)
    # 2nd column of solution is
    #   rho*E = rho*(Cv*T + u^2/2) = p/(gam-1) + rho*M^2*a^2/2
    #         = p/(gam-1) + gam*M^2*p/2 = (p/p0in)*p0in*(1/(gam-1) + gam*M^2/2)
    soln[:,2] = ps_by_p0in*p0in*(1./(gam-1) + (gam/2)*np.square(Ms))
    return soln
#enddef quasi_1D_steady_soln_in


def quasi_1D_steady_soln_residual(soln,prblm_setup,xarr,Dx=None):
    """
    Evaluates the residual field of the governing equations
    
    INPUTS:
    soln        : Solution in conserved variables from 'quasi_1D_steady_soln_in'
    prblm_setup : Problem setup data structure as output from call to function
                  'quasi_1D_steady_prblm_setup'
    xarr        : Array of grid coordinates along axis of nozzle
    Dx          : 1st-order finite difference operator
    
    OUTPUTS:
    Res : 2D array with x-grid points along the rows and the residuals of the
          mass, momentum & energy govening equations along the columns
    """
    # Evaluate the area ratio at the grid points
    Aarr = quasi_1D_steady_A_ovr_At(prblm_setup,xarr)
    # Extract conserved quantities as individual solution variables
    rho = soln[:,0];    rhou = soln[:,1];   rhoE = soln[:,2]
    # Derive primitive variables from the conserved variables in the solution
    _, u, p = quasi_1D_steady_soln_primitive(soln,prblm_setup)
    if Dx is None:  #Finite difference operator is not provided
        # Form the finite difference operator for 1st-order derivative of
        # 4th-order accuracy
        Dx = Findiff_Taylor_uniform(len(xarr),1,4)/(xarr[1]-xarr[0])
    # Compose the residual array
    Res = np.zeros_like(soln)  #Pre-allocate the residual array
    Res[:,0] = Dx.dot(rhou*Aarr)   #Residual of mass cons. eqn.
    Res[:,1] = Dx.dot((rhou*u+p)*Aarr)-p*Dx.dot(Aarr)  #Residual of mom. eqn.
    Res[:,2] = Dx.dot((rhoE + p)*u*Aarr)   #Residual of energy eqn.
    return Res
#enddef quasi_1D_steady_soln_residual


def quasi_1D_steady_soln_bcs(soln,prblm_setup):
    """
    Evaluates the solution of the quasi-1D steady Euler equations for the C-D
    nozzle problem corresponding to the boundary conditions (inlet density, and
    inlet & outlet pressure)
    
    INPUTS:
    soln        : Solution in conserved variables from 'quasi_1D_steady_soln_in'
    prblm_setup : Problem setup data structure as output from call to function
                  'quasi_1D_steady_prblm_setup'
    
    OUTPUTS:
    rhoin : Inlet density
    pin   : Inlet pressure
    pout  : Outlet pressure
    """
    # Derive primitive variables from the conserved variables in the solution
    rho, _, p = quasi_1D_steady_soln_primitive(soln,prblm_setup)
    # Return the values of inlet density and inlet & outlet pressure
    return rho[0], p[0], p[-1]
#enddef quasi_1D_steady_soln_bcs


def quasi_1D_steady_soln_primitive(soln,prblm_setup):
    """
    The solution is available in conserved variables; convert it to primitive
    variables
        
    INPUT:
    soln        : Solution in conserved variables from 'quasi_1D_steady_soln_in'
    prblm_setup : Problem setup data structure as output from call to function
                  'quasi_1D_steady_prblm_setup'

    OUTPUT:
    rho : Density, extracted from 0th column of 'soln' [1D array; length of soln]
    u   : Velocity
    p   : Pressure
    """
    # Parse the problem setup
    gam = prblm_setup[6]        #Sp. heat ratio
    # Extract 'rho'
    rho = soln[:,0]
    # Derive 'u' from 'rhou' as well as the 'rho' already extract above
    u = soln[:,1]/rho
    # Derive 'p' from 'rhoE' as well as the 'rho' and 'u' already obtained
    p = (soln[:,2] - 0.5*rho*np.square(u))*(gam-1)
    return rho, u, p
#enddef quasi_1D_steady_soln_primitive


def quasi_1D_steady_soln_plot(xarr,solns_in,add_in=None,add_title=None,\
        legends=None):
    """
    Plot solution of the quasi 1D steady Euler problem for C-D nozzles

    INPUTS:
    xarr      : Array of grid coordinates along axis of nozzle
    solns_in  : List of 2D arrays corresponding to solutions to plot
    add_in    : One optional additional field (same in no. as solns_in) to plot
    add_title : Optional title of subplot showing the additional field
    legends   : Optional ;egend strings corresponding to the solutions to plot

    OUTPUTS:
    None
    """
    if isinstance(solns_in,np.ndarray): #Numpy ndarray supplied
        assert solns_in.ndim in [2,3], 'solns_in must have 2 or 3 dimensions'
        # Form 'solns' as a list of solutions to plot, from given 'solns_in'
        if solns_in.ndim == 2:  #2D array, signifying a single solution to plot
            solns = [solns_in]  #Form a singleton list with it
        else:   #3D array given with different solutions to plot in the 3rd dim
            # Form list with entries as the 2D arrays corresponding to each
            # solution to plot
            solns = [solns_in[:,:,im] for im in range(solns_in.shape[2])]
    elif isinstance(solns_in,list):     #List supplied
        solns = solns_in    #Nothing to do; just copy
    vars = ['rho','rhou','rhoE']            #Variable names in solution
    colors = ['b','g','r','c','m','y','k']  #Various colors to use
    linetypes = ['-','--','-.',':']         #Various linetypes to use
    nSolns = len(solns)                     #No. of solution to plot
    plt.figure()
    for iv in range(len(vars)): #Each solution variable in its own subplot
        plt.subplot(2,2,iv+1)
        for iSoln in range(nSolns): #Go thru all solutions
            # We first vary the color for successive solutions; if we run out
            # of various colors, then we vary the linetype; this way we can
            # plot len(colors)*len(linetypes) solutions in distinct styles
            color = colors[iSoln % len(colors)] #Choose color; avoid repetition
            linetype = linetypes[iSoln // len(colors)]  #Choose linetype also
            plt.plot(xarr,solns[iSoln][:,iv],color+linetype) #Actual plot
        plt.title(vars[iv]) #Write the variable name as title of subplot
    plt.xlabel('x') #3rd subplot (bottom left) gets xlabel as it's in bottom row
    if add_in is not None:  #Additional field is to be plotted
        # Similar to obtain 'solns' from 'solns_in', obtain 'add' from 'add_in'
        if isinstance(add_in,np.ndarray):
            assert add_in.ndim in [2,3], 'add_in must have 2 or 3 dimensions'
            if add_in.ndim == 1:
                add = [add_in]
            else:
                add = [add_in[:,im] for im in range(add_in.shape[1])]
        elif isinstance(add_in,list):
            add = add_in
        assert len(add)==nSolns,'add_in should have as many solutions as solns'
        # Plot additional solution field just like the main solution, but in the
        # 4th subplot
        plt.subplot(2,2,4)
        for iSoln in range(nSolns):
            color = colors[iSoln % len(colors)]
            linetype = linetypes[iSoln // len(colors)]
            plt.plot(xarr,add[iSoln],color+linetype)
        plt.xlabel('x') #This one gets an xlabel as it's in bottom row
        if add_title is not None:   #Title is provided for additional field
            plt.title(add_title)
    if legends is not None: #Legend is given for solutions; add to last subplot
        plt.legend(legends)
#enddef quasi_1D_steady_soln_plot
