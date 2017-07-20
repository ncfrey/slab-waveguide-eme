#Author: Nathan Frey
#Date Created: 7.3.14

"""
Version History 1.0 7.17.14
Create GUI for numerical eigenmode expansion and optimization

8.5.14
GUI for calling eigenmode expansion calculations and simple optimization
"""

import Numerical_EME as NEME
from tkinter import *
from tkinter import ttk
import matplotlib.pyplot as plt
import numpy as np


#Nothing displayed yet, global reference to root window
plot_exists = False
global root



#Main Tk window


def ngui():

    #make gui window
    global root
    root = Tk()
    root.title("Numerical EME")

    #frame to put widgets in
    mainframe = ttk.Frame(root, padding="5 5 12 12")
    mainframe.grid(column=0, row=0, sticky=(N,W,E,S))
    mainframe.columnconfigure(0,weight=1)
    mainframe.rowconfigure(0,weight=1)

    #inputs (set global so other functions can grab them)
    global widths,nArray,mesh,z_len,lambda_0,T,Amp,nsVar,msVar,plVar
    widths = StringVar() #widths of each slab

    nArray = StringVar() #refractive indices of slabs

    mesh = StringVar() #mesh size
    z_len = StringVar() #z dimension of waveguide
    lambda_0 = StringVar() #wavelength of incoming light
    T = StringVar() #period of pert
    Amp = StringVar() #amplitude of pert

    #Option Menu inputs
    nsVar = StringVar(mainframe)
    nsVar.set("Number of slabs:")
    nsOptions = ["1","2","3","4","5","6","7","8","9"]
    numSlabs = OptionMenu(mainframe, nsVar, *nsOptions) #number of slabs in waveguide
    numSlabs.pack()
    numSlabs.grid(column=1,row=1,sticky=W)

    msVar = StringVar(mainframe)
    msVar.set("Mesh size:")
    msOptions = ["128","256","512","1024"]
    mesh = OptionMenu(mainframe,msVar,*msOptions) #Mesh size in x & z directions
    mesh.pack() 
    mesh.grid(column=1,row=2,sticky=W)

    plVar = StringVar(mainframe)
    plVar.set("Perturbed Layer:")
    plOptions = nsOptions                
    pLayer = OptionMenu(mainframe,plVar,*plOptions)
    pLayer.pack()
    pLayer.grid(column=1,row=3,sticky=W)

    #boxes for inputs
    widths_entry = ttk.Entry(mainframe,width=14,textvariable=widths)
    nArray_entry = ttk.Entry(mainframe,width=14,textvariable=nArray)
    z_len_entry = ttk.Entry(mainframe,width=7,textvariable=z_len)
    lambda_0_entry = ttk.Entry(mainframe,width=7,textvariable=lambda_0)
    T_entry = ttk.Entry(mainframe,width=7,textvariable=T)
    Amp_entry = ttk.Entry(mainframe,width=7,textvariable=Amp)

    #set default var values for standard OPV struct
    widths.set(".2, .1, .04, .14, .14, .1")
    nArray.set("1.5+7.8j, 2.0, 1.4+.1j, 1.7, 1.5, 1.0") #bottom to top
    z_len.set("10")
    lambda_0.set(".650")
    T.set("5")
    Amp.set(".02")
    msVar.set("512")
    nsVar.set("6")
    plVar.set("3")

    #Default values for single-slab glass waveguide (simplest case)
##    widths.set("1,1,1")
##    nArray.set("1.0,1.5,1.0")
##    z_len.set("100")
##    lambda_0.set(".650")
##    T.set("100*pi")
##    Amp.set(".25")
##    msVar.set("512")
##    nsVar.set("3")
##    plVar.set("2")
    

    #labels for inputs
    ttk.Label(mainframe,text="Widths").grid(column=2,row=1,sticky=(W,E))
    ttk.Label(mainframe,text="Refractive Indices").grid(column=2,row=2,sticky=(W,E))
    ttk.Label(mainframe,text="z length").grid(column=2,row=3,sticky=(W,E))
    ttk.Label(mainframe,text="lambda_0").grid(column=4,row=3,sticky=(W,E))
    ttk.Label(mainframe,text="Pert Period").grid(column=2,row=4,sticky=(W,E))
    ttk.Label(mainframe,text="Pert Amp").grid(column=4,row=4,sticky=(W,E))

    #buttons for calculations & help text
    ttk.Button(mainframe,text="Eigenmode Expansion",command=do_eme).grid(column=1,row=5,sticky=W)
    ttk.Button(mainframe,text="HELP",command=ngui_help).grid(column=5,row=5,sticky=W)
    ttk.Button(mainframe,text="Plot GM",command=plot_gm).grid(column=2,row=5,sticky=W)
    ttk.Button(mainframe,text="Plot Power", command=plot_pwr).grid(column=3,row=5,sticky=W)
    ttk.Button(mainframe,text="Show Structure",command=show_struct).grid(column=4,row=5,sticky=W)
    ttk.Button(mainframe,text="Quit",command=do_quit).grid(column=5,row=1,sticky=W)
    ttk.Button(mainframe,text="Optimize",command=do_opt).grid(column=1,row=4,sticky=W)

    #arrange input boxes
    widths_entry.grid(column=3,row=1,sticky=(W,E))
    nArray_entry.grid(column=3,row=2,sticky=(W,E))

    z_len_entry.grid(column=3,row=3,sticky=(W,E))
    lambda_0_entry.grid(column=5,row=3,sticky=(W,E))
    T_entry.grid(column=3,row=4,sticky=(W,E))
    Amp_entry.grid(column=5,row=4,sticky=(W,E))


    #padding/spacing around each input
    for child in mainframe.winfo_children():child.grid_configure(padx=5, pady=5)
    
    #Initialize GUI
    root.mainloop()

def do_eme():

    plt.close('all')
    widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val = get_values()
    callCode = "doeme"
    NEME.make_WG(callCode,widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val)
    
def plot_pwr():

    plt.close('all')
    widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val = get_values()
    callCode = "plotpwr"
    NEME.make_WG(callCode,widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val)


    plt.show()
    

def plot_gm():

    plt.close('all')
    widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val = get_values()
    callCode = "plotgm"
    NEME.make_WG(callCode,widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val)
    

    plt.show()
    
    
    

def show_struct():

    plt.close('all')
    #get user inputs
    widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val = get_values()

    #Send to Numerical_EME!!!
    callCode = "showstruct"
    fig = plt.figure()
    NEME.make_WG(callCode,widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val)

    plt.show()

#Run optimization
def do_opt():

    #grab user inputs as initial guesses
    widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val = get_values()

    callCode = "doopt"
    iter_num = 30
    fitness = 0
    R = 0

    width = widthsList[pLayer_val-2] #thickness of perturbed layer

    #Broad-band spectrum from 400 to 800 nm
    #Run loop for lambda_0 in bb to do broad band optimization
    bb = np.linspace(.400,.800,10)

    #initial mutation strengths
    o_width = width / 1000
    o_amp = Amp_val / 1000
    o_T = T_val / 1000

    #initial values
    w0 = width
    a0 = Amp_val
    t0 = T_val
    
    #best params
    _width = width
    _amp = Amp_val
    _T = T_val

    i = 0
    #Absorption over generations
    R_array = np.zeros((1))
    while i < iter_num:


        #R is overlap integral between incident E(x,0) and guided modes
        #A_tot contains absorption in each layer
        #Want to maximize this value for a given structure
        R, A_tot = NEME.make_WG(callCode,widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val)


        #Let's look at absorption fraction in absorber layer
        #as merit function - fifth layer in our OPV struct
        R = A_tot[4]
        R_array = np.insert(R_array,-1,R)
     
        
        
        #update fitness value if we have better absorption
        if R > fitness:
            fitness = R
            _width = width
            _amp = Amp_val
            _T = T_val

        #Revert to best parent
        if R < fitness:
            width = _width
            Amp_val = _amp
            T_val = _T

        

        #update values by mutation
        mutants = evolve(width,Amp_val,T_val,o_width,o_amp,o_T)
        o_width = mutants[0,0]
        del_w = mutants[0,1]
        o_amp = mutants[1,0]
        del_amp = mutants[1,1]
        o_T = mutants[2,0]
        del_T = mutants[2,1]
        

        width += del_w
        Amp_val += del_amp
        T_val += del_T
        widthsList[pLayer_val-2] = width

        
        #Check that values are within reasonable bounds
        if width < .8*w0 or width > 1.2*w0:
            width = _width
        if Amp_val < .5*a0 or Amp_val > 1.5*a0:
            Amp_val = _amp
        if T_val < .5*t0 or T_val > 1.5*t0:
            T_val = _T
        

        i += 1


    R_array = sorted(R_array)
    plt.plot(R_array)
    plt.title('Standard OPV Optimization')
    plt.xlabel('Number of generations')
    plt.ylabel('Guided mode overlap')
    plt.show()
    
    #Display optimized parameters
    print("Losses:",np.real(1-fitness),"Layer width:",_width,"Perturbation amplitude:",_amp,
          "Perturbation period:",_T)
    


#evolve parent structure
def evolve(width,Amp_val,T_val,o_width,o_amp,o_T):


    
    tau = 1 / np.sqrt(3) # 1 / sqrt of number of variable parameters degrees of freedom
    eps = np.random.normal(0,1,1) #a random, normally distributed variable

    #mutate generation only after first iteration
    if o_width != width / 100:
        
        o_width = o_width * np.exp(tau*eps) #log-normal (self-adaptive) update of mutation strength
        o_amp = o_amp * np.exp(tau*eps)
        o_T = o_T * np.exp(tau*eps)

##    del_w = np.random.normal(width,o_width,1)
##    del_amp = np.random.normal(Amp_val,o_amp,1)
##    del_T = np.random.normal(T_val,o_T,1)

    del_w = np.random.normal(0,o_width,1)
    del_amp = np.random.normal(0,o_amp,1)
    del_T = np.random.normal(0,o_T,1)

    #Let's try larger mutations
    del_w = 10*del_w
    del_amp = 10*del_amp
    del_T = 10*del_T

    

    print(del_w,del_amp,del_T)

    
    #array to hold mutated values
    mutants = np.zeros((3,2))
    mutants[0,0] = o_width
    mutants[0,1] = del_w
    mutants[1,0] = o_amp
    mutants[1,1] = del_amp
    mutants[2,0] = o_T
    mutants[2,1] = del_T
    
    

    return mutants
    
    
    

#This code block will collect user inputs for use in all other functions
def get_values():

    #Collect inputs and convert to appropriate data type
    numSlabs_val = nsVar.get()
    nSlabs = int(numSlabs_val)
    widths_val = widths.get()
    n_val = nArray.get()

    mesh_val = int(msVar.get())
    lambda_val = float(lambda_0.get())

    if plVar.get() != "Perturbed Layer:": #only get these values if user specifies pert
        pLayer_val = int(plVar.get())
        z_val = int(z_len.get())
        Amp_val = float(Amp.get())
        T_val = T.get()
        if T_val.endswith("pi"): #Handle n*pi period
            nCoeff = float(T_val.split("*")[0])
            T_val = nCoeff * np.pi
        elif T_val.endswith("pi") is False:
            T_val = float(T.get())

    #Period = (z_len*pi) / T
    ##T_val = (z_val*np.pi) / T_val
    
    #Parse widths and n's
    widthsList = widths_val.split(",")
    nList = n_val.split(",")

    #Allow only integer numSlabs number of widths and refractive indices
    widthsList = widthsList[0:nSlabs]
    nList = nList[0:nSlabs]
    
    #Convert strings to floats and complexes
    for i in range(0,nSlabs):
        widthsList[i] = float(widthsList[i])
        nList[i] = eval(nList[i])

    #Return 'em
    return widthsList,nList,mesh_val,lambda_val,z_val,Amp_val,T_val,pLayer_val
  
        







#Display help text for instructions on formatting of inputs, etc.
def ngui_help():

    title = "EME Help"
    message = ("Please input wavelength of incident light in microns (10^-6 m), "
               "e.g., 500 nm => .5. "
               "Likewise, input slab widths, length in z direction, "
               "perturbation amplitude, and perturbation period, in microns. "
               "Separate width and refractive index entries with commas. "
               " Contact: nfrey213@gmail.com with any questions.")

    messagebox.showinfo(title,message)

#Quit button
def do_quit():
    root.quit()
    root.destroy()


        
        

    




    
