from refl1d.names import *
from numpy import *


# "probe" refers to the measurement we will analyze, essentially telling Refl1D to 
# load a 4-column dataset with Q_Z, Reflectivity, dReflectivity, and dQ_Z. We can also
# add in some sample broadening and a variety of other parameters.
probe = list(range(2))
probe[0] = load4('YIG-Py_700mT.refl', sample_broadening=0.026)
probe[1] = load4('YIG-Py_4mT.refl', sample_broadening=0.026)
# It should be noted that sample broadening is an extremely important parameter for determining the
# shape of the critical edge, and it most often be determined through a bit of trial and error.  In
# theory, the broadening is the amount that the curvature of the film and/or substrate surface broadens
# the beam beyond the width one would expect based purely on the instrument resolution for a given
# slit configuration.
    
  
# Define Initial Beam Parameters
I0 = 1.00 # normalization factor
th = 0.000 # theta offset in degrees


#=========================================================================================================
# Define various parameters we are going to use later for better organization. 
#=========================================================================================================
# We will group these parameters by the material or layer we use them in. Generally "XX_rho" refers to the
# nuclear scattering length density (SLD) of material XX. Similarly, "XX_t" refers to layer thickness and
# "XX_sig" refers to interface roughness (both in Å). "XX_rhoM" will be magnetic SLD (To convert 
# magnetization to proper SLD units, take value in emu/cc * 0.00291). Of course, these names are arbitrary
# and you can use what's best for you.

#Ru
Ru_rho = 5.191
Ru_t = 30 
Ru_sig = 5


#Py (Ni80Fe20)
Py_rho = 9.15
Py_t = 200
Py_sig = 5
Py_rhoM = 2.46

#YIG
YIG_rho = 5.845
YIG_t = 400
YIG_sig = 5
YIG_rhoM = 0.4

#Pt
Pt_rho = 6.357
Pt_t = 100
Pt_sig = 4

# SiO2
SiO_t= 497.5
SiO_sig= 5
SiO_rho = 3.545

# Si
Si_sig = 5.97 
Si_rho = 2.069

# Remember, the above values are just our initial guesses - they don't have to be perfect. However, we
# recommend that you fix some values that you know, such as the nuclear SLD of the substrate.


#=========================================================================================================
# Define the actual materials which will make up our layers
#=========================================================================================================
Ru = SLD(name="Ru",rho = Ru_rho)
Py = SLD(name="Py",rho = Py_rho)
YIG = SLD(name="YIG",rho = YIG_rho)
Pt = SLD(name="Pt",rho = Pt_rho)
SiO= SLD(name="SiO",rho=SiO_rho)
Si= SLD(name="Si",rho=Si_rho)


#=========================================================================================================
# Build the sample for a variety of material layers
#=========================================================================================================
sample = list(range(2))
for i in range(len(sample)):
    sample[i] = (Si(100,Si_sig)
                |SiO(SiO_t,SiO_sig)
                |Pt(Pt_t,Pt_sig)
                |YIG(YIG_t,YIG_sig,Magnetism(YIG_rhoM,270, dead_above = 0, dead_below = 0))
                |Py(Py_t,Py_sig, Magnetism(Py_rhoM,270, dead_above = 0, dead_below = 0))
                |Ru(Ru_t,Ru_sig)
                |air)
			
# The sample is our stack of Si/SiO2/Pt/Y3Fe5O12/Ni80Fe20/Ru/air. Note that the format is (LayerA(LayerA_thickness,LayerA_interface)|LayerB(LayerB_thickness,LayerB_interface)|air
# If the layer is magnetic, you can instead use (LayerA(LayerA_thickness,LayerA_interface,Magnetism(rhoM = XXXX, SomeOtherParameter = XXX))|LayerB...|air)
# Please note that there are many magnetic parameters, including dead_above, dead_below, interface_above, interface_below, etc...

		
#=========================================================================================================
# Constraints 
#=========================================================================================================
# Here I like to put in any constraint that will help limit the model fitting. For example, we know that the neutron beam must have the
# same intensity, background, and angle offset for all scattering cross sections. The ".mm" and ".pp" refers to the "minus-minus" and
# "plus-plus" cross sections. Of course, the spin flip would then be ".pm" and ".mp". probe.mm.theta_offset = probe.pp.theta_offset.
# We can also accomplish this be including a line which says "probe.shared_beam()". 
for i in range(len(probe)):
    probe[i].mm.intensity = probe[i].pp.intensity
    probe[i].mm.theta_offset = probe[i].pp.theta_offset
    probe[i].mm.background = probe[i].pp.background
    probe[i].mm.sample_broadening = probe[i].pp.sample_broadening

# Here we are going to link the thicknesses and interface roughnesses of each layer to be the same across the two measurements. The only
# differences are going to be in the magnetism of the samples. We don't need to explicitly link the nuclear SLDs because they are defined
# by the "material", so if we use the same materials in each sample construction it'll be automatically linked by default.
for i in range(0,5):
    sample[1][i].thickness = sample[0][i].thickness
    sample[1][i].interface = sample[0][i].interface

# Here we initialize the starting values for the incoming neutron beam, background, and sample misalignment.
for i in range(len(probe)):
    probe[i].pp.theta_offset.value = th
    probe[i].pp.intensity.value = I0

#=========================================================================================================
# Fit parameters 
#=========================================================================================================
# Now we will begin selecting which parameters ought to be fitted, and what range they should vary over. To do this
# we just say sample[XX][LayerXX].parameterXX.range(bottomRange,topRange). We could also say ".pm(XX)" or ".pmp(XX)" instead
# of ".range(XX,YY)" if we want to have the parameter vary by ±XX or ±XX%, respectively.

# Here we fit the beam intensity, in case our normalization is a bit off (often not necessary) and the theta offset, in case our
# alignment is a tiny bit off. One can also bit a constant background, but this is generally discouraged since we have explicitly
# measured and subtracted off the background.
for i in range(len(probe)):
    probe[i].pp.intensity.range(0.9,1.1)
    probe[i].pp.theta_offset.range(-0.005,0.005)
    probe[i].pp.sample_broadening.pmp(30)	

# Here we have the option to fit the thicknesses
if 1:
    sample[0][Ru].thickness.range(15,45)
    sample[0][Py].thickness.range(150,250)
    sample[0][YIG].thickness.range(300,450)
    sample[0][Pt].thickness.range(75,125)
    sample[0][SiO].thickness.range(300,600)

# Here we have the option to fit the interface roughnesses It's important to note that the roughness parameter 
# captures  long-range roughness contributions such as thickness variation and conformal roughness as well as 
# chemical intermixing.
if 1:
    sample[0][Ru].interface.range(1,20)
    sample[0][Py].interface.range(1,20)
    sample[0][YIG].interface.range(1,30)
    sample[0][Pt].interface.range(1,40)
    sample[0][SiO].interface.range(2,30)
    sample[0][Si].interface.range(1,60)

# Here we have the option to fit the nuclear SLDs. Note that we are not fitting the Si SLD since I know the crystal is
# high quality so that I can input theoretical density and composition values into the SLD calculator at
# https://www.ncnr.nist.gov/resources/activation/. So we don't fit things we know if it's not necessary to describe the
# data. On the other hand, we don't know the actual composition or density of the YIG film, so we may fit the SLD.
if 1:
    Ru.rho.range(5,6)
    Py.rho.pmp(8)
    YIG.rho.pmp(20) #.pmp(20) means fit range of 20% around initial value
    Pt.rho.pmp(5)
    SiO.rho.range(2,4.1)

# Here we have to option to fit various magnetic properties
if 1:
    sample[0][Py].magnetism.rhoM.pmp(20)
    sample[1][Py].magnetism.rhoM.pmp(20)
    sample[0][Py].magnetism.dead_above.range(0,20)
    sample[1][Py].magnetism.dead_above.range(0,20)
    sample[0][YIG].magnetism.rhoM.range(0,2)
    sample[1][YIG].magnetism.rhoM.range(-2,2)
    sample[0][YIG].magnetism.dead_above.range(0,40)
    sample[1][YIG].magnetism.dead_above.range(0,40)
# Magnetic dead layers, which have the same nuclear SLD as rest of the layer but a magnetization of 0
# can be added at the top and bottom of each layer


#=========================================================================================================
# Fitting Problem - technical definition 
#=========================================================================================================
zed = 0.5 # microslabbing bin size, in Å
alpha = 0.0 # integration factor - leave this at zero
step = False 
#  This "step" or "step_interfaces" setting is really important. Effectively, there are two ways to calculate the reflectivity. 
# One way treats the entire system as a series of microslabe ("step_interfaces = True") while the other uses something
# called the Nevot-Croce approximation ("step_interfaces = False"). The difference between the two is beyond the scope of these
# notes, but generally you should know that Nevot-Croce is much, much faster but fails in some very specific cases where
# the sample contains layers for which the roughness on one side is very different from the other AND the roughness is comparable
# to or larger than the thickness. Be careful here , and always use "step_interfaces=False" as a sanity check.

# Define a single model object for each measurement we want Refl1D to fit. We give it the sample structure (sample), the dataset to fit (probe), etc.
models = []
model0 = Experiment(sample=sample[0], probe=probe[0], dz=zed,dA=alpha,step_interfaces=step)
models.append(model0)
model1 = Experiment(sample=sample[1], probe=probe[1], dz=zed,dA=alpha,step_interfaces=step)
models.append(model1)
#Now tell Refl1D to actually fit the data:
problem = MultiFitProblem(models)
problem.name = "YIG-Py"