# === Import section ===
import sys
import numpy as np
from molgroups import mol
from molgroups import components as cmp
from molgroups import lipids
from refl1d.names import Parameter, SLD, Slab, Experiment, FitProblem, load4
from refl1d.probe import QProbe
from refl1d.flayer import FunctionalProfile
from periodictable.fasta import Sequence, H2O_SLD, D2O_SLD

# === Constant definition section ===
# Canvas
DIMENSION = 300
STEPSIZE = 0.5

# Hermite Spline
CONTROLPOINTS = 7
SPACING = 15.0

# SLDS
NSLDH2O = H2O_SLD
NSLDD2O = D2O_SLD

def bilayer(z, sigma, bulknsld, global_rough, rho_substrate, l_lipid1, l_lipid2, l_tether,  nf_tether, mult_tether, vf_bilayer):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld = bulknsld * 1e-6
    rho_substrate = rho_substrate * 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=global_rough, rho_substrate=rho_substrate,
              l_lipid1=l_lipid1, l_lipid2=l_lipid2, l_tether=l_tether, nf_tether=nf_tether, mult_tether=mult_tether,
              vf_bilayer=vf_bilayer, radius_defect=1e8)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = blm.fnWriteProfile(z)
    normarea = blm.normarea

    # for statistical analysis of molgroups
    problem.moldat, problem.results = write_groups([blm], ['bilayer'])
    
    # Return nSLD profile in Refl1D units
    return apply_bulknsld(z, bulknsld, normarea, area, nsl) * 1e6

def write_groups(groups, labels):
    """Return dictionaries with combined output of fnWriteGroup2Dict and fnWriteResults2Dict
    
        Inputs:
        groups: list of Molgroups objects to process
        labels: list (same length as groups) of labels"""
    
    moldict = {}
    resdict = {}
    for lbl, gp in zip(labels, groups):
        moldict = {**moldict, **gp.fnWriteGroup2Dict({}, lbl, np.arange(DIMENSION) * STEPSIZE)}
        resdict = {**resdict, **gp.fnWriteResults2Dict({}, lbl)}
        
    return moldict, resdict

def apply_bulknsld(z, bulknsld, normarea, area, nsl):
    """Given area and nSL profiles, fill in the remaining volume with bulk material"""
    
    # Fill in the remaining volume with buffer of appropriate nSLD
    nsld = nsl / (normarea * np.gradient(z)) + (1.0 - area / normarea) * bulknsld

    # Return nSLD profile in Refl1D units
    return nsld

def make_samples(func, substrate, contrasts, **kwargs):
    """Create samples from combining a substrate stack with a molgroups layer
    
        Inputs:
        func: function used to define FunctionalProfile object. Must have form func(z, bulknsld, *args)
        substrate: Refl1D Stack or Layer object representing the substrate
        contrasts: list of buffer materials, e.g. [d2o, h2o]. One sample will be created for each contrast
        **kwargs: keyword arguments. Must have one keyword argument for each arg in func(..., *args), but
                  not one for bulknsld"""
    samples = []

    for contrast in contrasts:
        mollayer = FunctionalProfile(DIMENSION * STEPSIZE, 0, profile=func, bulknsld=contrast.rho, **kwargs)
        layer_contrast = Slab(material=contrast, thickness=0.0000, interface=5.0000)
        samples.append(substrate | mollayer | layer_contrast)

    return samples

# Define bilayer parameters
vf_bilayer = Parameter(name='volume fraction bilayer', value=0.9).range(0.0, 1.0)
l_lipid1 = Parameter(name='inner acyl chain thickness', value=10.0).range(8, 30)
l_lipid2 = Parameter(name='outer acyl chain thickness', value=10.0).range(8, 18)
sigma = Parameter(name='bilayer roughness', value=5).range(2, 9)
global_rough = Parameter(name ='substrate roughness', value=5).range(2, 9)
d_oxide = Parameter(name='silicon oxide layer thickness', value=10).range(5, 30)
d_Cr =  Parameter(name='chromium layer thickness', value=40).range(10, 150)
d_gold =  Parameter(name='gold layer thickness', value=140).range(100, 200) #thickness of gold
rough_cr_au =  Parameter(name='gold chromium roughness', value=10).range(2, 24.0) # roughness of Cr/Au interface
l_tether =  Parameter(name='tether thickness', value=10).range(3, 50) #distance from substrate to inner headgroup/acyl chain interface
nf_tether = Parameter(name='number fraction tether', value=0.5).range(0.1, 0.7)
mult_tether = Parameter(name='multiplicity tether', value=1.0).range(0.1, 10)

### Define bilayer object
#POPS = cmp.Lipid(name='POPS', headgroup=cmp.ps, tails=[cmp.palmitoyl, cmp.palmitoyl], methyls=[cmp.methyl])
#POPE = cmp.Lipid(name='POPE', headgroup=cmp.pe, tails=[cmp.palmitoyl, cmp.oleoyl], methyls=[cmp.methyl])
#DOPIP2 = cmp.Lipid(name='DOPIP2', headgroup=cmp.pip2, tails=[cmp.oleoyl, cmp.oleoyl], methyls=[cmp.methyl])

blm = mol.tBLM(tether=lipids.WC14, filler=cmp.bme, lipids=[lipids.DOPC],
                lipid_nf=[1.0])

## === Stack ===
##
## First, we create a 'material' for each bulk layer, which has an real and imaginary
## scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='d2o', rho=6.3000, irho=0.0000)
h2o = SLD(name='h2o', rho=-0.56, irho=0.0000)

# nSLD parameters
d2o.rho.range(5.3000, 6.5000)
h2o.rho.range(-0.6, 2.6)

## Substrate materials
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)

siox = SLD(name='siox', rho=4.1000, irho=0.0000)
siox.rho.range(2.9000, 5.1000)

cr = SLD(name='chromium', rho=2.7, irho=0.0)
cr.rho.range(2.2000, 4.0000)

gold = SLD(name='gold', rho=4.4, irho=0.0) #iro is the absorption of neutrons, should be 0
gold.rho.range(4.2000, 4.8000)

## Then bulk layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:

layer_silicon = Slab(material=silicon, thickness=0.0000, interface=global_rough)
layer_siox = Slab(material=siox, thickness=d_oxide, interface=global_rough)
layer_cr = Slab(material=cr, thickness=d_Cr, interface=rough_cr_au)
layer_gold = Slab(material=gold, thickness=d_gold - (blm.substrate.z + 0.5 * blm.substrate.length), interface=0.0000)

## Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.

substrate = layer_silicon  | layer_siox | layer_cr | layer_gold

sample, sampleh = make_samples(bilayer, substrate, [d2o, h2o], sigma=sigma, 
                             global_rough=global_rough, rho_substrate=gold.rho,
                             l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
                             nf_tether=nf_tether, mult_tether=mult_tether,
                             vf_bilayer=vf_bilayer)


## === Data files ===
datapath = './dat/'

# substrate only
# probe = load4(datapath + 'cr001_d2o.refl', back_reflectivity=True, name='D2O buffer')
# probeh = load4(datapath + 'cr002_h2o.refl', back_reflectivity=True, name='H2O buffer')
probe = load4(datapath + 'cr001_d2o.ort', back_reflectivity=True, name='D2O buffer')
probeh = load4(datapath + 'cr002_h2o.ort', back_reflectivity=True, name='H2O buffer')

# Background parameter
# probe.background.value = 0.0000
probe.background.range(-1e-7, 1e-5)
probeh.background.range(-1e-7, 1e-5)

probe.intensity.range(0.90, 1.15)
probeh.intensity = probe.intensity
probe.sample_broadening.range(-0.003, 0.01)
probe.theta_offset.range(-0.01, 0.02)

# Define critical edge oversampling for samples that require it (typically D2O only)
probe.critical_edge(substrate=silicon, surface=d2o)

## === Problem definition ===
## a model object consists of a sample and a probe,

## step = True corresponds to a calculation of the reflectivity from an actual profile
## with microslabbed interfaces.  When step = False, the Nevot-Croce
## approximation is used to account for roughness.  This approximation speeds up
## the calculation tremendously, and is reasonably accuarate as long as the
## roughness is much less than the layer thickness
step = False

model = Experiment(sample=sample, probe=probe, dz=STEPSIZE, step_interfaces = step)
modelh = Experiment(sample=sampleh, probe=probeh, dz=STEPSIZE, step_interfaces = step)

#problem = FitProblem([model0, model0h, model, model_prot, modelh_prot])
problem = FitProblem([model, modelh])
#problem = FitProblem([model, modelh])

def custom_plot():

    import plotly.graph_objs as go
    from refl1d.webview.server.colors import COLORS

    moldat = problem.moldat

    def hex_to_rgb(hex_string):
        r_hex = hex_string[1:3]
        g_hex = hex_string[3:5]
        b_hex = hex_string[5:7]
        return int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)

    n_lipids = 1
    group_names = {'gold substrate': ['bilayer.substrate'],
               #'silicon oxide': ['bilayer.siox'],
               'bME': ['bilayer.bME'],
               'tether': ['bilayer.tether_bme', 'bilayer.tether_free', 'bilayer.tether_hg'],
               'tether acyl chains': ['bilayer.tether_methylene', 'bilayer.tether_methyl'],
               'inner headgroups': [f'bilayer.headgroup1_{i}' for i in range(1, n_lipids + 1)],
               'inner acyl chains': [f'bilayer.methylene1_{i}' for i in range(1, n_lipids + 1)] + [f'bilayer.methyl1_{i}' for i in range(1, n_lipids + 1)],
               'outer acyl chains': [f'bilayer.methylene2_{i}' for i in range(1, n_lipids + 1)] + [f'bilayer.methyl2_{i}' for i in range(1, n_lipids + 1)],
               'outer headgroups': [f'bilayer.headgroup2_{i}' for i in range(1, n_lipids + 1)],
              }
    
    normarea = moldat['bilayer.normarea']['area']

    fig = go.Figure()
    traces = []
    MOD_COLORS = COLORS[1:]
    color_idx = 1
    sumarea = 0
    for lbl, item in group_names.items():
        area = 0
        for gp in item:
            if gp in moldat.keys():
                zaxis = moldat[gp]['zaxis']
                area += moldat[gp]['area']
            else:
                print(f'Warning: {gp} not found')

        color = MOD_COLORS[color_idx % len(MOD_COLORS)]
        plotly_color = ','.join(map(str, hex_to_rgb(color)))
        traces.append(go.Scatter(x=zaxis,
                                 y=area / normarea,
                                 mode='lines',
                                 name=lbl,
                                 line=dict(color=color)))
        traces.append(go.Scatter(x=zaxis,
                                 y=area / normarea,
                                 mode='lines',
                                 line=dict(width=0),
                                 fill='tozeroy',
                                 fillcolor=f'rgba({plotly_color},0.3)',
                                 showlegend=False
                                 ))
        color_idx += 1
        sumarea += area

    color = COLORS[0]
    plotly_color = ','.join(map(str, hex_to_rgb(color)))
    
    traces.append(go.Scatter(x=zaxis,
                                y=sumarea / normarea,
                                mode='lines',
                                name='buffer',
                                line=dict(color=color)))
    traces.append(go.Scatter(x=zaxis,
                                y=sumarea / normarea,
                                mode='lines',
                                line=dict(width=0),
                                fill='tonexty',
                                fillcolor=f'rgba({plotly_color},0.3)',
                                showlegend=False
                                ))    
    traces.append(go.Scatter(x=zaxis,
                                y=[1.0] * len(zaxis),
                                mode='lines',
                                line=dict(color=color, width=0),
                                showlegend=False))

    
    fig.add_traces(traces[::-1])

    fig.update_layout(
        title='Component Volume Occupancy',
        template = 'plotly_white',
        xaxis_title=dict(text='z (Ang)'),
        yaxis_title=dict(text='volume occupancy')
    )

    return fig

setattr(problem, 'custom_plot', custom_plot)

if __name__ == '__main__':

    problem.chisq_str()

    moldat = problem.moldat

    fig = custom_plot(moldat)
    fig.show()