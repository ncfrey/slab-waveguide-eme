#Author: Nathan Frey
#Date Created: 6.19.14

"""
Version History
1.0 6.26.14
Non-perturbative treatment of multilayer slab waveguide w/ single frequency
incident light.

2.0 7.18.14
Incorporate small dielectric perturbation (surface corrugation)
Build in GUI
Calculate E(x,0) from incoming plane wave

8.5.14
Complete eigenmode expansion treatment of complex-indexed waveguides

TO DO:
*Include TM modes...
*Check TMM, initial and final A,B coeffs, proper shape of fields
*To get from Aj to Aj+1, must use the Trans Mat between these two layers, not
 full Trans Mat
"""

from numpy import *
from numpy.linalg import eig, solve, norm, inv
from matplotlib.pyplot import *
from numpy.ma import conjugate
from scipy.integrate import simps, trapz, odeint
from scipy import signal
import cmath
import math



    
#Initialize waveguide struct, calculate modes, incident beam, do expansion
def make_WG(callCode,wid, nArray, nx, lambda_0, z_len, S_0, p_period, pLayer):
    
    """
    Inputs:
    callCode tells which routine called the function and what routines
    should be run
    wid is array with thickness of each layer
    nx is number of mesh points
    nArray is refractive index in each slab
    lambda_0 is wavelength of incoming light

    z_len is length in z dimension (propagation direction)
    del_e is dielectric perturbation
    p_period is period of perturbation
    pLayer is layer where perturbation occurs

    Outputs:
    Supported guided mode eigenvalues (propagation constants) & eigenmodes
    Eigenmode decomposition of incident E field in terms of gm & all modes (leaky + guided)
    Losses due to leaky modes (given by overlap integral)
    Coupling coefficients for the guided modes    
    """

    nz = nx #z mesh size
    numSlabs = len(wid)
    d = cumsum(wid) #thickness of entire structure        
    wg = linspace(0, d[numSlabs - 1], nx) #mesh of structure
    c = 3e14 #in microns
    k = (2*pi) / lambda_0 #wavenumber
    omega_0 = k*c

    #construct z profile
    z = linspace(0,z_len,nz)

    #construct perturbation, del_e
    T = (z_len*pi) / p_period #period == (z_len*pi) / T
    phi = S_0 #y offset 
    del_e = S_0 * signal.square(z*T) + phi


    dx = d[numSlabs - 1] / nx #dx in the waveguide between mesh points

    #construct unperturbed eps_0 dielectric profile in x direction
    eps_0 = ones((nx),'complex')
    for i in range(0,nx):
        pos = i*dx
        if pos < d[1]:
            eps_0[i] = nArray[0]**2
        for j in range(1,numSlabs): #figure out which slab we're in
            if pos >= d[j - 1] and pos < d[j]:
                eps_0[i] = nArray[j]**2

    

 

    
    eps = zeros((nx,nz),'complex')  #construct permittivity profile in x direction
    eps_plt = zeros((nx,nz)) #real part of eps, for plotting purposes
    for i in range(0,nx):
        pos = i*dx
        if pos < d[1]:
            eps[i,:] = nArray[0]**2
            eps_plt[i,:] = real(nArray[0])**2
        for j in range(1,numSlabs): #figure out which slab we're in
            if pos >= d[j - 1] and pos < d[j]:
                eps[i,:] = nArray[j]**2
                eps_plt[i,:] = real(nArray[j])**2
                


    #characterize perturbation
    pmax = amax(del_e)


    #construct permittivity in z direction
    p_layer = pLayer
    for i in range(0,nx):
        pos = i*dx
        if pos >= d[p_layer - 2] and pos < d[p_layer - 2] + pmax: #pos is in the layer
            for m in range(0,nz):
                if del_e[m] == pmax:
                    if p_layer - 1 > 0:
                        eps[i,m] = real(nArray[p_layer - 2])**2 #magnitude of dielectric perturbation
                        eps_plt[i,m] = real(nArray[p_layer - 2])**2
        

    
    #For perturbed bottom or top layer
    #Make periodic gaps filled with air (plasmon resonance)
    for i in range(0,nx):
        pos = i*dx
        if p_layer - 1 <= 0: #looking at bottom layer
            if pos >= 0 and pos < d[p_layer - 1]:
                for m in range(0,nz):
                    if del_e[m] == pmax:
                        eps[i,m] = 1.0 #Air
                        eps_plt[i,m] = 1.0
        if p_layer == numSlabs: #looking at top layer
            if pos >= d[p_layer - 2] and pos < d[-1]:
                for m in range(0,nz):
                    if del_e[m] == pmax:
                        eps[i,m] = 1.0 #Air
                        eps_plt[i,m] = 1.0

    
            


    #Calculate wavelength in perturbed layer, diffraction angles, & intensities
    lambda_p, theta, I = lambda_calc(nArray,lambda_0,z_len,pLayer,p_period)

    #Get incident electric field for incoming plane wave @ lambda_0
    #Using unperturbed matrix method get
    E_struct,M_TE,M_TM,AB_0 = plane_wave(wid, nx, nArray, k, theta, I, pLayer)


    #Weird periodic layer electric field function - useless
    ##pl_field(pLayer,wid,nx,k,eps,p_period,z_len,E_struct,M_TE,M_TM,AB_0,lambda_p,theta)
                    
    ##ib = incident_beam(wg, k) #define incident electric field profile as gaussian beam

    ib = E_struct #incident beam E(x,z=0)
##    plot(real(ib),wg)
##    title('Incident Plane Wave Electric Field')
##    xlabel('Field distribution')
##    ylabel('Structure thickness')
##    show()
    
    #Figure out which module called me up w/ callCode
    #plot image of dielectric profile in structure
    if callCode == "showstruct":
        X, Z = meshgrid(z,wg)
        eps_img = pcolormesh(X,Z,real(eps_plt),cmap='Greens') #eps_plt is only real part of RIP
        colorbar(eps_img,orientation='horizontal')
        title('Dielectric profile')
    

    #Plot guided modes within the structure
    if callCode == "plotgm":
        ev, em = calc_modes(wg, k, eps_0, wid)
        gm, ee, c = EME(wg, ib, ev, em, eps_0,callCode)
    
    #Plot power propagation for simple structures
    if callCode == "plotpwr":
        E_xz = CN(wg, eps_0, k, ib, z_len, nz)

    #Return values for various calculations, including coupling coefficients
    if callCode == "doeme":
        ev, em = calc_modes(wg, k, eps_0, wid)
        gm, ee, c, c_frac, R = EME(wg, ib, ev, em, eps_0,"doeme")
        A_tot = calc_absorp(gm, wid)
        betas, A_integrated = calc_Axz(nArray, p_layer, gm, ee, k, del_e, p_period, z, c)

        betas = real(betas) #discard 0 complex parts
        A_integrated = real(A_integrated)

        data = vstack([betas,c_frac])
        data = vstack([data,A_integrated])
        headers = ["Propagation constants","Coupling coefficients","A(x,z) integrated"]
        for i in range(0,len(headers)):
            print(headers[i],": ",data[i])
        print("Absorption in layers: ", A_tot)
        print("Losses: ", real(1-R))

    #Perform simple evolution strategy optimization
    if callCode == "doopt":
        ev, em = calc_modes(wg, k, eps_0, wid)
        R, gm = EME(wg, ib, ev, em, eps_0,"doopt")

        #We want to look at absorption over ALL modes
        A_tot = calc_absorp(em, wid)
        ## A_tot = calc_absorp(gm,wid)
        return R, A_tot

    #For testing w/o using the gui
    if callCode == "slabwg":
        ##X, Z = meshgrid(z,wg)
        ##eps_img = pcolormesh(X,Z,real(eps_plt),cmap='Greens')
        ##show()
        ##E_xz = CN(wg, eps_0, k, ib, z_len, nz)
        ev, em = calc_modes(wg, k, eps_0, wid)
                
        gm, ee, c, c_frac, R = EME(wg, ib, ev, em, eps_0,"doeme")
        beta_tm0, TM0 = TM0_calc(wid, nx, nArray, k)
        
        plot(TM0,wg)
        show()
        A_tot = calc_absorp(TM0,wid)
        
    
        ##A_tot = calc_absorp(gm, wid)
        ##print(A_tot)
        
        ##plot(wg,gm)
        ##show()
    


        
        ##print(A_tot)
        

        
        
    #Mode calculations, unperturbed solutions
    ##ev, em = calc_modes(wg, k, eps_0, wid) #Get unperturbed eigenmodes, propagation constants
    ##gm, ee, c = EME(wg, ib, ev, em, eps_0,callCode) #do eigenmode expansion, return unperturbed guided modes
    ##A_tot = calc_absorp(gm, wid) #array with absorption % in each layer

    
    # A(x,z) via Yeh method
    ##calc_Axz(gm, ee, k, del_e, p_period, z, c)
    
    #Coeffs for entire (x,z) struct
    ##calc_cz(wg, eps, ib, k, z_len, nz)
    
    #Crank-Nicholson scheme for coupled coefficients
    #returns electric field distribution, E(x,z) thru entire structure
    ##E_xz = CN(wg, eps_0, k, ib, z_len, nz)


    
#Define params for example waveguide struct of interest
def slab_wg():

    #x,y dependent params
    ##wid = [4.0, 1.5, 1.0, 1.5, 1.0, 1.5, 1.5, 2.0] #normalized to wavelength
    nx = 512 #mesh size
    ##nArray = [1.49, 1.52+1j, 1.49, 1.52, 1.49, 1.52, 1.49, 1.0] #refractive indices in each layer
    lambda_0 = .50001 #normalized to 632 nm

    #Simple slab waveguide
    ##wid = [1.0,2.0,1.0]
    ##nArray = [1.0,1.5,1.0]


    #OPV structs from Yutong's IEEE paper
    #a)
    wid =[.2, .04, .01, .04, .14, .5, .1]
    nArray = [.8+6.1j, 3.5, 1.9+.4j, 3.5, 1.9, 1.5, 1.0]
        
    #z params
    z_len = 100
    nz = nx
    z = linspace(0,z_len,nz)

    #perturbation
    p_layer = int(floor(len(wid)/2)) #layer to apply perturbation to
    T = (z_len / 10) * pi #period == (z_len*pi) / T
    S_0 = .2*wid[p_layer] #Amplitude = 2*S_0
    phi = S_0 #y offset
    
    del_e = S_0 * signal.square(z*T) + phi
    p_period = (z_len*pi) / T



    #OPV params
    ##z_len = 10
    ##lambda_0 = .650
    ##p_period = 5
    ##S_0 = .02
    
    callCode = "slabwg"
    make_WG(callCode,wid, nArray, nx, lambda_0, z_len, S_0, p_period,p_layer)
    

#calculate all eigenvalues (effective permittivities)
#and corresponding eigenmodes
def calc_modes(wg, k, eps, wid):

    nx = len(wg)
    dx = wg[1] - wg[0]
    c = 3e14 #microns/s


    eps = conjugate(eps)
    eps = transpose(eps)


    #calc matrix of hermitian wave operator L
    Ljj = array(((-2 * ones(nx) / (dx**2 * k**2))))
    maindiag = array((Ljj + eps))
    
    offdiag = array(((1*ones(nx - 1) / (dx**2 * k**2))))

    L = diag(offdiag,-1)+ diag(maindiag,0) + diag(offdiag,1)
    
    


    #ev are eigenvalues, effective permitivites (beta / k)^2
    #em are corresponding eigenmodes
    ev, em = eig(L)

    



    return ev, em

#Mode condition matrix method calculation TM0 mode
def TM0_calc(wid, nx, nArray, k):

    lambda_0 = (2*pi) / k #wavelength of incoming light
    nLayers = len(wid)
    d = cumsum(wid)
    beta_tm0 = 0

    AB = zeros((nLayers,2),'complex') #Coeffs for TM field
    ns = nArray[0] #substrate index of refraction
    n0 = nArray[nLayers - 1] #incident medium index (air)
    
    #Storage needed for TMM
    k_x = zeros((nLayers),'complex')
    cos_theta = zeros((nLayers),'complex')
    phi = zeros((nLayers),'complex')
    D = zeros((nLayers,4),'complex')
    P = zeros((nLayers,4),'complex')

    #Window for effective epsilons for guided modes
    epseff_0 = real(conjugate(n0)*n0)
    epseff_s = real(conjugate(ns)*ns)
    dEff = .01
    

    for eps_eff in arange(epseff_0 + dEff,epseff_s - dEff,(epseff_s - epseff_0 - 2*dEff)/100000):
        for l in range(0,nLayers):

            #Terms for TMM
            beta = k*sqrt(eps_eff)
            k_x[l] = cmath.sqrt( k**2 * conjugate(nArray[l])*nArray[l] - beta**2)
            cos_theta[l] = k_x[l] / (k*nArray[l])
            nx_l = (wid[l] / d[nLayers-1]) * nx
            phi[l] = k_x[l] * nx_l

            #Dynamical matrix for TM waves
            D[l,0] = cos_theta[l]
            D[l,1] = cos_theta[l]
            D[l,2] = nArray[l]
            D[l,3] = -nArray[l]

            #Fill in propagation matrix
            P[l,0] = exp(1j*phi[l])
            P[l,1] = 0
            P[l,2] = 0
            P[l,3] = exp(-1j*phi[l])

        transfer_mat = zeros((2,2),'complex')
        D_l = zeros((2,2),'complex')
        P_l = zeros((2,2),'complex')
        
        for l in range (1,nLayers-2):
            D_l[0,0] = D[l,0]
            D_l[0,1] = D[l,1]
            D_l[1,0] = D[l,2]
            D_l[1,1] = D[l,3]
            
            
            P_l[0,0] = P[l,0]
            P_l[0,1] = P[l,1]
            P_l[1,0] = P[l,2]
            P_l[1,1] = P[l,3]

            if l == 1:
                transfer_mat = D_l.dot(P_l).dot(inv(D_l))

            if l > 1:
                A = D_l.dot(P_l).dot(inv(D_l))
                transfer_mat = dot(transfer_mat,A)

            tm = transfer_mat

##            D_l[0,0] = D[nLayers-1,0]
##            D_l[0,1] = D[nLayers-1,1]
##            D_l[1,0] = D[nLayers-1,2]
##            D_l[1,1] = D[nLayers-1,3]
##
##           
##            tm = dot(inv(D_l),tm)

##            #substrate
##            D_l[0,0] = D[0,0]
##            D_l[0,1] = D[0,1]
##            D_l[1,0] = D[0,2]
##            D_l[1,1] = D[0,3]
##
##            #completed transfer matrix
##            tm = dot(tm,D_l)
            
            AB[0,0] = 0
            AB[0,1] = 1
            ##AB[l,0] = dot(tm,AB[l-1,0])[0,0]
            ##AB[l,1] = dot(tm,AB[l-1,1])[1,0]
            AB[l,:] = dot(tm,AB[l-1,:])
            
        #incident medium
        D_l[0,0] = D[nLayers-1,0]
        D_l[0,1] = D[nLayers-1,1]
        D_l[1,0] = D[nLayers-1,2]
        D_l[1,1] = D[nLayers-1,3]

       
        transfer_mat = dot(inv(D_l),transfer_mat)

        #substrate
        D_l[0,0] = D[0,0]
        D_l[0,1] = D[0,1]
        D_l[1,0] = D[0,2]
        D_l[1,1] = D[0,3]

        #completed transfer matrix
        transfer_mat = dot(transfer_mat,D_l)
        M_TM = transfer_mat

        AB[nLayers-1,:] = dot(M_TM,AB[nLayers-2,:])
        
        ##AB[nLayers-1,0] = dot(M_TM, AB[nLayers-2,0])[0,0]
        ##AB[nLayers-1,1] = dot(M_TM, AB[nLayers-2,0])[1,0]


        #Mode condition for TM modes is that M[0,0] == 0
        #Give a little wiggle room for computational error
        if abs(real(M_TM[0,0])) < .001:
            M_TM[0,0] = 0 #correct computational error
            beta_tm0 = beta


            #Calc TM field for this beta
            E_m = zeros((),'complex')
         

            #Get A and B constant coeffs for each layer
            for m in range(0,nLayers):
                #Initial conditions
                if m == 0:
                    AB[0,0] = 0 #A_0 for TM confined mode that vanish at infinity
                    AB[0,1] = 1 #B_0
                    nx_m = ((d[0]) / d[nLayers-1]) * nx

                elif m > 0:
                    ##AB[m,:] = dot(M_TM,AB[m-1,:]) #A_m & B_m

                    nx_m = ((d[nLayers-m]-d[nLayers-m-1])/d[nLayers-1]) * nx #num pts in layer

                
                #Deal with float nx_m
                nx_m = round(nx_m)
                    
                x_m = linspace(0,wid[m],nx_m) #x space in incident layer
                E_m = AB[m,0] * exp(-1j*k_x[m]*x_m) + AB[m,1] * exp(1j*k_x[m]*x_m)

                if m == 0:
                    E_struct = E_m #Incident electric field @ E(x,z=0)
                elif m > 0:
                    E_struct = append(E_struct,E_m)


            
            #Trim or pad edges
            if len(real(E_struct)) > nx:
                overflow = len(E_struct) - nx
                for i in range(0,overflow):
                    E_struct = delete(E_struct,i)
            if len(real(E_struct)) < nx:
                while len(E_struct) < nx:
                    E_struct = insert(E_struct,-1,E_struct[len(E_struct)-1])

            
            norm = amax(E_struct) #normalize to 1
            TM0 = E_struct / norm
            return beta_tm0, TM0
            break


    #No mode found
    if beta_tm0 == 0:
        return 0, 0

#define E(x,0) based on incoming plane wave using TMM
def plane_wave(wid, nx, nArray, k, theta, I, pLayer):

    
    lambda_0 = (2*pi) / k #wavelength of incoming light
    nLayers = len(wid)
    d = cumsum(wid)

    ns = nArray[0] #substrate index of refraction
    n0 = nArray[nLayers - 1] #incident medium index (air)

    
    

    #Calculate terms needed for TMM
    k_x = zeros((nLayers),'complex')
    cos_theta = zeros((nLayers),'complex')
    phi = zeros((nLayers),'complex')
    D = zeros((nLayers,4),'complex')
    P = zeros((nLayers,4),'complex')

    #Total E field over all betas
    E_tot = zeros((nx,len(theta)),'complex')

    #A, B coeffs for electric field at each layer
    AB = zeros((nLayers,4),'complex')
    AB[0,0] = 1 #A_0 for TE
    AB[0,1] = 0 #B_0
    AB[0,2] = 0 #For TM
    AB[0,3] = 1
           


    #Calc E at every beta value, integrate over intensities of diffraction angles
    for nn in range(0,len(theta)):
        #Calculate transfer matrix for TE & TM waves
        for i in range(0,2):        
            for l in range(0,nLayers):

                #Beta is 0 until light hits perturbed layer (zero refraction, no z propagation)
                #perturbation -> non-zero beta

                #propagation constant, wavenumber in layer, cosine of ray angle
                #num pts in the layer, phase
                beta = cmath.sqrt(k**2 * conjugate(nArray[l])*nArray[l] - \
                                  k**2 * conjugate(nArray[l])*nArray[l]*cos(theta[nn]))
                k_x[l] = cmath.sqrt( k**2 * conjugate(nArray[l])*nArray[l] - beta**2)
                cos_theta[l] = real(k_x[l] / (k*nArray[l]))
                nx_l = (wid[l] / d[nLayers-1]) * nx
                phi[l] = k_x[l] * nx_l

                

                
                #Fill in dynamical matrix for TE waves
                if i == 0:
                    D[l,0] = 1
                    D[l,1] = 1
                    D[l,2] = nArray[l]*cos_theta[l]
                    D[l,3] = -nArray[l]*cos_theta[l]

                #Dynamical matrix for TM waves
                if i == 1:
                    D[l,0] = cos_theta[l]
                    D[l,1] = cos_theta[l]
                    D[l,2] = nArray[l]
                    D[l,3] = -nArray[l]

                #Fill in propagation matrix
                P[l,0] = exp(1j*phi[l])
                P[l,1] = 0
                P[l,2] = 0
                P[l,3] = exp(-1j*phi[l])

            cos_s = real(cos_theta[0]) #for use in transmittance calculation
            cos_0 = real(cos_theta[nLayers - 1])


            transfer_mat = zeros((2,2),'complex')
            D_l = zeros((2,2),'complex')
            P_l = zeros((2,2),'complex')
            
            for l in range (1,nLayers-2):
                D_l[0,0] = D[l,0]
                D_l[0,1] = D[l,1]
                D_l[1,0] = D[l,2]
                D_l[1,1] = D[l,3]

                P_l[0,0] = P[l,0]
                P_l[0,1] = P[l,1]
                P_l[1,0] = P[l,2]
                P_l[1,1] = P[l,3]

                if l == 1:                    
                    transfer_mat = D_l.dot(P_l).dot(inv(D_l))

                    if i == 0:
                        M_TE = transfer_mat
                        AB[l,0:2] = dot(M_TE,AB[l-1,0:2]) #A_m & B_m for TE

                    if i == 1:
                        M_TM = transfer_mat
                        AB[l,2:4] = dot(M_TM,AB[l-1,2:4]) #A_m & B_m for TM
                        
                if l > 1:
                    A = D_l.dot(P_l).dot(inv(D_l))
                    transfer_mat = dot(transfer_mat,A)

                    if i == 0:
                        M_TE = transfer_mat
                        AB[l,0:2] = dot(M_TE,AB[l-1,0:2]) #A_m & B_m for TE

                    if i == 1:
                        M_TM = transfer_mat
                        AB[l,2:4] = dot(M_TM,AB[l-1,2:4]) #A_m & B_m for TM
                    

                #Must apply TMM at this step to get coeffs for each layer
                    
            #incident medium
            D_l[0,0] = D[nLayers-1,0]
            D_l[0,1] = D[nLayers-1,1]
            D_l[1,0] = D[nLayers-1,2]
            D_l[1,1] = D[nLayers-1,3]

            transfer_mat = dot(inv(D_l),transfer_mat)

            #substrate
            D_l[0,0] = D[0,0]
            D_l[0,1] = D[0,1]
            D_l[1,0] = D[0,2]
            D_l[1,1] = D[0,3]

            #completed transfer matrix
            transfer_mat = dot(transfer_mat,D_l)

            #Store
            if i == 0:
                M_TE = transfer_mat
                AB[-1,0:2] = dot(M_TE,AB[nLayers-3,0:2])

            if i == 1:
                M_TM = transfer_mat
                AB[-1,2:4] = dot(M_TM,AB[nLayers-3,2:4])


        
        ##for M in (M_TE,M_TM):
            ##R = real((conjugate(M[1,0])*M[1,0]) / (conjugate(M[0,0])*M[0,0])) #Reflectance

            ##if R > 1:
                ##R = 1.0 #correct for slight numerical error
                   
            ##if cos_0 > 0:
                ##T = real((ns*cos_s)/(n0*cos_0) * (1/(conjugate(M[0,0])*M[0,0]))) #Transmittance           

                ##print("T", T, "R", R, M)

        #Set up incident plane wave electric field

        d = cumsum(wid) #thickness of structure
        E_m = zeros((),'complex')
     

        #Get A and B constant coeffs for each layer
        for m in range(0,nLayers):
##            #Initial conditions
##            if m == 0:
##                AB[0,0] = 1 #A_0 for TE
##                AB[0,1] = 0 #B_0
##                AB[0,2] = 0 #For TM
##                AB[0,3] = 1
##                nx_m = ((d[0]) / d[nLayers-1]) * nx
##
##            elif m > 0:
##                AB[m,0:2] = dot(M_TE,AB[m-1,0:2]) #A_m & B_m for TE
##                AB[m,2:4] = dot(M_TM,AB[m-1,2:4]) #A_m & B_m for TM
                
##                AB[m,0] = dot(M_TE,AB[m-1,0])[0,0] #A_m
##                AB[m,1] = dot(M_TE,AB[m-1,1])[1,0] #B_m
##                AB[m,2] = dot(M_TM,AB[m-1,2])[0,0]
##                AB[m,3] = dot(M_TM,AB[m-1,3])[1,0]
            
##                nx_m = ((d[nLayers-m]-d[nLayers-m-1])/d[nLayers-1]) * nx #num pts in layer


        #Get number of points in each layer
            if m == 0:
                nx_m = ((d[0]) / d[nLayers-1]) * nx
            elif m > 0:
                nx_m = ((d[nLayers-m]-d[nLayers-m-1])/d[nLayers-1]) * nx
        
            #Deal with float nx_m
            nx_m = round(nx_m)
                
            x_m = linspace(0,wid[m],nx_m) #x space in incident layer
            E_m = AB[m,0] * exp(-1j*k_x[m]*x_m) + AB[m,1] * exp(1j*k_x[m]*x_m)

            if m == 0:
                E_struct = E_m #Incident electric field @ E(x,z=0)
            elif m > 0:
                E_struct = append(E_struct,E_m)

        
        #Trim or pad edges
        if len(real(E_struct)) > nx:
            overflow = len(E_struct) - nx
            for i in range(0,overflow):
                E_struct = delete(E_struct,i)
        if len(real(E_struct)) < nx:
            while len(E_struct) < nx:
                E_struct = insert(E_struct,-1,E_struct[len(E_struct)-1])
                                     

        #E_struct is unperturbed solution, now include diffraction angles (non-zero beta)
        #due to periodic perturbation
        
        
        E_tot[:,nn] = (E_struct) #throw out time-dependent phase (imaginary part)


    #E_struct is weighted average of E @ each beta from diffraction orders
    #weighted by intensity of diffraction order
    nn = 0
    E_struct = zeros((nx),'complex')
    
    for nn in range(0,len(theta)):
        E_struct += I[nn]*E_tot[:,nn]



    return E_struct, M_TE, M_TM, AB




#Calculate wavelength and diffracted orders of incoming light in perturbed layer
def lambda_calc(nArray,lambda_0, z_len, pLayer,p_period):

    pLayer -= 1 #Subtract one so pLayer matches with 0 indexing of nArray
    numSlabs = len(nArray)
    l_m = lambda_0
    
    for m in range(numSlabs-1,pLayer-1,-1):
        lambda_p = lambda_0 / real(nArray[m]) #real part of n measures refraction
    theta = zeros(()) #0th order diffraction = 0
   
    m = 1
    #Loop until all diffraction orders are calculated
    while m >= 0:
        if (m*lambda_p/p_period) > 1: #break from loop once all orders are calculated
            break
        theta_m = arcsin((m*lambda_p/p_period)) #diffraction angles b/c of grating
        theta = append(theta,theta_m)
        m += 1

    I = zeros((len(theta))) #Intensity of each diffraction order
    m = 0
    N = z_len / p_period #number of slits
    delta = 2*pi*p_period / lambda_p #period dependent term

    #From diffraction grating intensity eq
    for m in range(0,len(theta)):
        if theta[m] == 0:
            I[m] = 0
        elif theta[m] > 0:
            I[m] = sin(N*delta*sin(theta[m]))**2 / (delta*sin(theta[m])/2)**2 

    norm = cumsum(I)[len(I)-1]
    I /= norm #Normalized to represent % intensities of each order
    
    
    
    return lambda_p, theta, I
    

#Incident electric Field E(x,z=0) in perturbed layer
def pl_field(p_layer,wid,nx,k,eps,p_period,z_len,E_struct,M_TE,M_TM,AB_0,lambda_p,theta):

    d = cumsum(wid) #total struct thickness
    dx = d[len(wid)-1] / nx #spacing per point
    nz = nx
    dz = z_len / nz
    pmax = 0
    eps_sum = 0
    point_count = 1
    
    
    for i in range(0,nx):
        pos = i*dx
        if pos >= d[p_layer - 2] and pos < d[p_layer - 1]: #we're in perturbed layer
            _pmax = real(cmath.sqrt(amax(eps[i,:])))
            if _pmax > pmax:
                pmax = _pmax #keep track of max eps value in the layer
            for m in range(0,nz):
                pos_z = m*dz
                eps_sum += eps[i,m]
                point_count += 1

    n0 = cmath.sqrt(eps_sum) / point_count #average index in the layer
    n1 = pmax - n0

    #terms to simplify coeff solutions
    k0 = k*n0
    K = (2*pi) / p_period
    del_k = 2*k0 - K
    kappa = k / (2*n1)
    s = cmath.sqrt( kappa**2 - (del_k/2)**2 )



    #Find unperturbed A,B constant coeffs in above layer
    if p_layer > 1:
        A_0 = AB[p_layer-1,0]
        B_0 = AB[p_layer-1,0]
    elif p_layer == 1:
        A_0 = 1
        B_0 = 1
    
        

    #Solutions for A(x), B(x) coefficients - from parabolic approximation
    x = linspace(0,d[len(wid)-1],nx)

    c1 = 0
    c2 = 0
    
    #c1 and c2 are constant coeffs that come from BC w/ other layers via TMM
    A = (c1*cosh(s*x) + c2*sinh(s*x))*exp(1j*(del_k/2)*x)
    B = exp(-1j*del_k*x/2)*(1j/kappa)*(sinh(s*x)*(1j*del_k*c2/2 + c1*s) \
                                       + cosh(s*x)*(1j*del_k*c1/2 + c2*s))

    #Electric field E(x,0) in the layer
    E = A*exp(-1j*k0*x) + B*exp(1j*k0*x)
    

    
        



#define an incident beam profile
def incident_beam(wg, k):
    nx = len(wg)
    dx = wg[1] - wg[0]
    


    f = zeros(nx)
    ##pos = 0.6
    ##wid = .5
    pos = wg[nx/2]
    wid = 10*dx
    
    E_0 = 1
    for m in range(0, nx):
        x = wg[m]
        #f[m] = E_0*cos(k*x) #can be modified to create different beam profile
        f[m] = E_0*exp(-((x-pos)/wid)**2) #gaussian beam profile
    f = conjugate(f)
    f = transpose(f)
    
  
    return f

#eigenmode expansion in terms of all eigenmodes (leaky + guided)
#and only guided; gives losses
def EME(wg, ib, ev, em, eps,callCode):

    cc = conjugate(em)
    cc = transpose(cc)

    
    c = dot(cc,ib) #coefficients for expansion
    f_all = dot(em,c) #expansion of incident E field in terms of all modes
    ##sigma_c = cumsum(sqrt(conjugate(c)*c))
    sigma_c = cumsum(c)

    #show the expansion in terms of all calculated modes
    ##plot(wg, f_all)
    ##show()

    
    
    nx = len(eps)
    ep_substrate = eps[0]
    ep_cover = eps[nx - 1]
    ep_thresh = amax([ep_substrate, ep_cover]) #guided modes must be greater than this threshold


    #find eigenvectors that correspond to guided mode eigenvalues
    gm = em[:, ev > ep_thresh]
    ee = ev[ev > ep_thresh]

    
    

    #Normalize to 1
    for i in range(0,gm.shape[1]):
        norm = amax(gm[:,i])
        gm[:,i] /= norm
  

    #calculate coefficients (c) for decomposition in terms of guided modes
    cc = conjugate(gm)
    cc = transpose(cc)
    c = dot(cc,ib)

    

    #gives percentage of coupling into each of the gm
    c_frac = real((conjugate(c)*c) / (conjugate((sigma_c[-1]))*sigma_c[-1]))
    

    f_gm = dot(gm, c) #expansion in terms of only guided modes
    norm = amax(f_gm)
    f_gm /= norm
    
    #show the expansion and/or the guided modes IFF any exist
    ##plot(wg, f_gm, wg, f_all)
    ##show()
    if gm.shape[0] > 0 and callCode == "plotgm":
        gmplot = plot(real(gm),wg)
        xlabel('Normalized field distribution')
        ylabel('Structure thickness')
        title('Guided modes')

##        fgmplot = plot(ib,wg,f_gm,wg)
##        ylabel('Structure thickness')
##        xlabel('Incident field distribution')
##        title('Guided Mode Expansion')

    #Normalize all vectors to avoid overflow error
    norm = amax(f_all)
    f_all /= norm
    norm = amax(ib)
    ib /= norm
    norm = amax(f_gm)
    f_gm /= norm


    #Overlap between all modes and incident electric field
    ##R_all = real(dot(conjugate(dot(conjugate(f_all),ib)),dot(conjugate(f_all),ib)) / \
    ##                           ( dot(conjugate(f_all),f_all) * dot(conjugate(ib),ib) ))

    R_all = real(dot(conjugate(dot(conjugate(f_all),ib)),dot(conjugate(f_all),ib)) / \
                               ( dot(conjugate(f_all),f_all) * dot(conjugate(ib),ib) ))

    R = real(dot(conjugate(dot(conjugate(f_gm),ib)),dot(conjugate(f_gm),ib)) / \
                                ( dot(conjugate(f_gm),f_gm) * dot(conjugate(ib),ib) ))

    #Is this necesary for complex guided modes?
    gm = sqrt(conjugate(gm)*gm)
    
    # 1 - R represents losses
    if callCode == "doopt":
        return R, gm

    if callCode == "plotgm" or callCode == "slabwg":
        return gm, ee, c_frac

    if callCode == "doeme":
        return gm, ee, c, c_frac, R

def calc_cz(wg, eps, ib, k, z_len, nz):

    #unperturbed eps
    eps = eps[:,1]
    ev, em = calc_modes(wg, k, eps)


    num = len(ev)
    pos_ev = zeros((num),'complex')
    count = 0

    #pick out only positive effective permittivities, i.e. forward propagating modes
    for m in range(0,num):
        if ev[m] > 0 :      
            pos_ev[count] = ev[m]
            count += 1
    pos_ev = trim_zeros(pos_ev)

    #forward propagating modes
    pos_m = em[:,ev > 0]

    #positive, real-valued propagation constants
    betas = k * sqrt(pos_ev)

    pos_m = conjugate(pos_m)
    pos_m = transpose(pos_m)

    z = linspace(0,z_len,nz)
    
    num = len(pos_ev)
    c = dot(pos_m,ib)
    cz = zeros((num,nz))
    for m in range(0, num):
        cz[m,:] = real(c[m]*exp(1j*betas[m]*z))

    pos_m = conjugate(pos_m)
    pos_m = transpose(pos_m)

    #pick out gm from pos_m

    nx = len(eps)
    gm_markers = zeros((nx))
    ep_substrate = eps[0]
    ep_cover = eps[nx - 1]
    ep_thresh = amax([ep_substrate, ep_cover]) #guided modes must be greater than this threshold


    #find eigenvectors that correspond to guided mode eigenvalues
    count = 0
    for m in range(0,num):
        if pos_ev[m] > ep_thresh:
            gm_markers[count] = m
            count += 1
    gm_markers = trim_zeros(gm_markers)
    gm_markers = gm_markers.astype(int)

    #A(x,z) coupling coeffs for guided modes
    gm_xz_coeffs = cz[gm_markers,:]
    gm_xz_coeffs = transpose(gm_xz_coeffs)
    

    
    ib_z = dot(pos_m, cz)
    X, Z = meshgrid(z,wg)
    #visualize propagation of power through waveguide
    contour(X,Z,abs(ib_z)**2,20)
    show()
    

    
    
#Crank-Nicholson scheme for calculating coupled coeffs A(x,z)
def CN(wg, eps, k, ib, z_len, nz):

    z = linspace(0, z_len, nz)
    dz = z[1] - z[0]
    n0 = cmath.sqrt(eps[0]) #refractive index of layer
    nx = len(wg)
    dx = wg[1] - wg[0]

    #Generate Propagation operator matrix (P)
    eps = conjugate(eps)
    eps = transpose(eps)
    
    maindiag = array( (-2*ones(nx)/dx**2/(2*k*n0)),'complex') #((For a single slab/substrate))
    ##maindiag = array(( (-2*ones(nx) / dx**2 + k**2 * (eps - n0**2)) / (2*k*n0) ),'complex')
    offdiag = array((ones(nx - 1) / dx**2 / (2 * k * n0)),'complex')
    P = array((diag(offdiag,-1) + diag(maindiag,0) + diag(offdiag,1)),'complex')


    #Forward & backward propagation
    Forw = identity((nx),'complex') + 1j * dz*P / 2
    Back = identity((nx),'complex') - 1j * dz*P / 2

 
    ib_z = zeros((nx,nz),'complex')
    thresh = 1e-5 #cutoff for electric field

    ib_z[:,0] = ib #incident field at E(x,z=0)

    #Solve for ib_z at each z value, modify according to Hadley BC
    #that is, transparent boundary for evanescent fields
    for m in range(1,nz):
        FF = Forw
        BB = Back

        
        #update fields at lower boundary
        if (abs(ib[0]) > thresh):

            k0 = 1j / dx * log(ib[1]/ib[0])
           
            if real(k0) <= 0:
                k0 = 1j * imag(k0)
            TBC = exp(1j*k0*dx) / dx**2 * 1j * dz / 2 / (2*k*n0)
            FF[0,0] = FF[0,0] + TBC
            BB[0,0] = BB[0,0] - TBC

        #and on upper boundary
        if (abs(ib[nx-1]) > thresh):
            k0 = -1j / dx * log(ib[nx-1] / ib[nx-2])
            if real(k0) <= 0:
                k0 = 1j * imag(k0)
            TBC = exp(1j*k0*dx) / dx**2 * 1j * dz / 2 / (2*k*n0)
            FF[nx-1,nx-1] = FF[nx-1,nx-1] + TBC
            BB[nx-1,nx-1] = BB[nx-1,nx-1] - TBC
        

        #solve for incident field at z with the transparent BC
        ib = solve(BB,dot(FF,ib))
        ib_z[:,m] = real(ib)
    
    X, Z = meshgrid(z,wg)
    #visualize propagation of power through waveguide
    pwrplt = contour(X,Z,abs(ib_z)**2,32)


    return ib_z

def calc_absorp(gm, wid):

   
    #Make everything positive
    gm = sqrt(conjugate(gm)*gm)

    nModes = len(gm.shape)
    nLayers = len(wid)
    
    if nModes > 1: #more than one guided mode (need to count columns)
        nModes = gm.shape[1]

    if nModes == 1:
        dy = len(gm[:])
    elif nModes > 1:
        dy = len(gm[:,0]) #number of points gm is evaluated at
    axis = cumsum(wid)
    axis = insert(axis,0,0)

    
    A = zeros((nLayers,nModes),'complex') #array holds absorption values in each layer for each mode
    A_s = zeros((nLayers,nModes),'complex')

    #Integrate guided modes to find absorption in each layer
    for m in range(0,nModes):
        for i in range(0,nLayers):
            nx = (axis[i+1] - axis[i]) / axis[-1] * dy #num of pts evaluated at within layer
            

            x0 = (axis[i] - axis[0]) / axis[-1] * dy #pt that layer starts at

            #Approximate by cutting off at float pt values
            nx = round(nx)
            x0 = round(x0)
            xf = x0 + nx #pt that layer ends at

##            xaxis = linspace(axis[i],axis[i+1],nx)
            xaxis = linspace(x0,xf,nx)
            
            ##a_t = trapz(gm[:,m],xaxis)
            if nModes == 1:
                a_s = simps(gm[x0:xf],xaxis) #use simpson's rule to approx integral
            elif nModes > 1:
                a_s = simps(gm[x0:xf,m],xaxis)
            ##A[i,m] = a_t #Trapezoidal & Simpson's rule agree
            A_s[i,m] = a_s

  
    A_tot = zeros((nLayers),'complex') #total absorption in a layer across all GMs
    for i in range(0,nLayers):
        a_tot = A_s[i,:].sum()
        A_tot[i] = a_tot
        
    #normalize layer absorption to 1
    norm = cumsum(A_tot)[nLayers - 1]
    A_tot = A_tot / norm #gives % absorption for GMs in each layer

##    A_gm = A_s #each column reps % of GM absorbed in a layer
##    for m in range(0,nModes):
##        norm = cumsum(A_gm[:,m])[nLayers - 1]
##        A_gm[:,m] = A_gm[:,m] / norm

    #A_tot returns absorption of layers in the order incident medium -> substrate (top -> bottom)
    return real(A_tot) #neglect very small imaginary part



def calc_Axz(nArray, p_layer, gm, ee, k, del_e, T, z, c):

    nModes = len(ee)
    betas = zeros((nModes),'complex')
    for i in range(0,nModes):
        betas[i] = real(k*cmath.sqrt(ee[i])) #propagation constants - discard 0 imaginary part

    if p_layer - 2 > 0:
        amp = real(nArray[p_layer - 2])**2 #magnitude of dielectric perturbation
    elif p_layer - 2 < 0:
        amp = 1.0
        
    c0 = 3*10**8 #speed of light in vacuum
    omega = k*c0 #freq of incoming light


    #need to do nCr (nModes Choose 2) number of calculations?
    #No - calculate A's in pairs

    A = zeros((len(z),nModes),'complex') #hold all A(z)
    A_int = zeros((nModes),'complex') #hold each |A(z)|**2

   
    for m in range(0,nModes-1):
        n = m+1
        
        gm_m = gm[:,m]
        gm_n = gm[:,n]
        c_m = c[m] #unperturbed coupling coeffs to be used as initial conditions
        c_n = c[n]

        gm_m = conjugate(gm_m)
        gm_m = transpose(gm_m)
        
        ovlap_mn = real(dot(gm_m,gm_n)) #overlap integral between modes - neglect small imaginary part
        del_beta = real(betas[m] - betas[n] - (2*pi / T)) #mode matching condition
        kappa_mn = (omega/4)*amp*ovlap_mn #represents coupling between modes
      

        #Soln for Am, An diff eqs comes from Mathematica
        #Need two coeffs, c1, c2, from initial conditions
        alpha_1 = (1j / (2*kappa_mn)) * (1j * del_beta - cmath.sqrt(-4*kappa_mn**2 - del_beta**2))
        alpha_2 = (1j / (2*kappa_mn)) * (1j * del_beta + cmath.sqrt(-4*kappa_mn**2 - del_beta**2))
        
        coeffMat = array(([1,1],[alpha_1, alpha_2]),'complex')
        solnMat = array(([c_m, c_n]),'complex')

        #Solve system of linear eqs for constant coeffs for A_m(z), A_n(z)
        new_c = solve(coeffMat, solnMat)
        c1 = (new_c[0]) #Real or complex?
        c2 = (new_c[1])

        A_m = array((len(z)),'complex')
        A_n = array((len(z)),'complex')

    

        #These expressions come from soln to {dAm/dz ~ An ; dAn/dz ~ Am}
        A_m = c1*exp( (-1/2)*z*(cmath.sqrt(-4*kappa_mn**2 - del_beta**2) - 1j*del_beta)) \
              + c2 * exp( (1/2)*z*(cmath.sqrt(-4*kappa_mn**2 - del_beta**2) + 1j*del_beta))

        A_n = (1j / (2*kappa_mn)) * ( c1*(1j*del_beta - cmath.sqrt(-4*kappa_mn**2 - del_beta**2)) \
                                      * exp( (-1/2)*z*(cmath.sqrt(-4*kappa_mn**2 - del_beta**2) + 1j*del_beta) ) \
                                      + c2*(cmath.sqrt(-4*kappa_mn**2 - del_beta**2) + 1j*del_beta) \
                                      * exp((1/2)*z*(cmath.sqrt(-4*kappa_mn**2 - del_beta**2)-1j*del_beta)) )

        A_m = (betas[m] / abs(betas[m])) * A_m #sign indicates propagation direction
        A_n = (betas[n] / abs(betas[n])) * A_n
        
        ccm = transpose(A_m)
        ccm = conjugate(ccm)
        ccn = transpose(A_n)
        ccn = conjugate(ccn)
        
        ##A[:,m] = real(A_m) #Am(z) values are complex...
        ##A[:,n] = real(A_n)

        A_int[m] = real(dot(ccm,A_m)) # |A(z)|**2
        A_int[n] = real(dot(ccn,A_n))

        ##print(real(A_int))
    norm = cumsum(real(A_int))[nModes-1]
    A_int = real(A_int / norm) #normalize coeffs to sum to 1
        
    return betas, A_int
    

    

        
            

            
            
            
            

            

    

        
        

            
        
    

    
    
    
    
    
    
        
    
    

    

    

    

    
    
    
