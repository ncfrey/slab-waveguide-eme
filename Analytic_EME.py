#Author: Nathan Frey
#Date Created: 6.24.14

"""
Version History:
1.0 8.5.14
Simple script for running root-finding based analytic calculations
of guided modes.
"""

from numpy import *
from scipy.optimize import fsolve, brentq
from matplotlib.pyplot import plot, show
from scipy import odeint

#make a simple slab waveguide

def solve_gm():
    """
    Defines simple slab waveguide, one layer, with incoming light at wavelength
    lambda_0, and solves transcendental eqs for propagation constants
    """

    n1 = 1.5
    n2 = 2.5
    lambda_0 = .5
    k = (2*pi) / lambda_0

    b_low = n1*k
    b_up = n2*k

    nb = 14 #mesh partition
    db = (b_up - b_low) / nb #mesh size

  
    x0 = zeros(nb)
    count = 0
    #b_low = int(floor(b_low))
    #b_up = int(ceil(b_up))
    
    #construct mesh
    for i in arange(b_low, b_up,db):
        x0[count] = i
        count += 1

    num = len(x0)
    sol_TE = zeros(nb)
    count = 0

    #check for sign change -> zero exists between mesh points
    for i in range(0,num - 1):
        if sign( f_TE(x0[i]) ) != sign (f_TE(x0[i+1])):
            sol_TE[count] = brentq(f_TE, x0[i],x0[i+1])
            count += 1

    sol_TE = trim_zeros(sol_TE)


    sol_TM = zeros(nb)
    count = 0
    for i in range(0,num - 1):
        if sign( f_TM(x0[i]) ) != sign (f_TM(x0[i+1])):
            sol_TE[count] = brentq(f_TM, x0[i],x0[i+1])
            count += 1

    sol_TM = trim_zeros(sol_TE)
    sol_GM = sol_TE / k #Guided mode betas

    
    return sol_GM
  

def f_TE(x):

    #define params
    n1 = 1.5 #air outside layers
    n2 = 2.5 #glass slab
    lambda_0 = .5 #incident light wavelength normalized to 500 nm
    k = (2*pi) / lambda_0 #wavenumber
    d = lambda_0 #slab thickness to 100 nm

    b_low = n1*k #bounds for guided modes
    b_up = n2*k
    
    #solve transcendental eqs for guided modes
    h = (b_up**2 - x**2)**(1/2)
    q = (x**2 - b_low**2)**(1/2)



    return tan(h*d) - (2*h*q) / (h**2 - q**2)

#Define incident electric field, E(x,0)
def incident_beam(d, k)

    dx = 512 #mesh size
    wg = linspace(0,d, dx)
    nx = len(wg)
    pos = 0
    wid = 1
    E_0 = 1

    ib = zeros(nx)
    for m in range(0,nx):
        x = wg[m]
        ib[m] = E_0*exp(-((x-pos)/wid)**2)

    return ib

def deriv(E,x, k, n2, betas):

    eps = sqrt(n2)
    num = len(betas)
    for m in range(0,num):
        beta_m = betas[m]
        return array(E[1],(beta_m**2 - k**2 * eps)*E[0])

def solve_e(deriv, wg):

    E_init = array([0,1])
    E = odeint(deriv,E_init, wg)
    

def f_TM(x):

    #define params
    n1 = 1.5 #air outside layers
    n2 = 2.5 #glass slab
    lambda_0 = .5 #incident light wavelength normalized to 500 nm
    k = (2*pi) / lambda_0 #wavenumber
    d = lambda_0 #slab thickness to 100 nm

    b_low = n1*k #bounds for guided modes
    b_up = n2*k

    h = (b_up**2 - x**2)**(1/2)
    qq = (n2**2 / n1**2) * (x**2 - b_low**2)**(1/2)

    return tan(h*d) - (2*h*qq) / (h**2 - qq**2)
    
