
import numpy as np
from molgroups import mol
from molgroups import components as cmp
from molgroups import lipids
from refl1d.names import load4, Parameter, SLD, Slab, Experiment, FitProblem
from refl1d.flayer import FunctionalProfile

### Define molgroups space.
DIMENSION=400       # Number of steps

# Length of steps. Also sets calculation resolution, and determines speed of calculation
STEPSIZE=0.5

## === Film structure definition section ===

### Bilayer profile definition function

def bilayer(z, sigma, bulknsld, global_rough, rho_substrate,nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld = bulknsld * 1e-6
    rho_substrate = rho_substrate * 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=sigma, rho_substrate=rho_substrate,
              nf_tether=nf_tether, mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
              vf_bilayer=vf_bilayer, radius_defect=1e8)

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = blm.fnWriteProfile(z)
    normarea = blm.normarea
    
    # for statistical analysis of molgroups
    problem.moldat, problem.results = write_groups([blm], ['bilayer'])
    
    # Return nSLD profile in Refl1D units
    return apply_bulknsld(z, bulknsld, normarea, area, nsl) * 1e6

def bilayer_prot(z, sigma, bulknsld, global_rough, rho_substrate,nf_tether, mult_tether, l_tether, l_lipid1, l_lipid2, vf_bilayer, gamma, beta, protz, protnf, protexchratio, radius_defects, l_mr, frac_pos_mr):
    """ Generic tethered bilayer """

    # Scale all SLDs from Refl1D units (1e-6 Ang^-2) to molgroups units (Ang^-2)
    bulknsld = bulknsld * 1e-6
    rho_substrate = rho_substrate * 1e-6

    blm.fnSet(sigma=sigma, bulknsld=bulknsld, global_rough=sigma, rho_substrate=rho_substrate,
              nf_tether=nf_tether, mult_tether=mult_tether, l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
              vf_bilayer=vf_bilayer, radius_defect=radius_defects)
    prot.protexchratio = protexchratio
    mr.fnSet(length=l_mr, frac_position=frac_pos_mr)
    protm.fnSet(gamma, beta, protz, sigma, protnf, bulknsld)    
    blmprot.fnAdjustBLMs()

    # Calculate scattering properties of volume occupied by bilayer
    normarea, area, nsl = blmprot.fnWriteProfile(z)
    normarea = blm.normarea

    # export objects for post analysis, needs to be from this function
    problem.moldat, problem.results = write_groups([blm, protm], ['bilayer','protein'])

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
vf_bilayer = Parameter(name='volume fraction bilayer', value=0.95).range(0.0, 1.0)
l_lipid1 = Parameter(name='inner acyl chain thickness', value=12.0).range(8, 30)
l_lipid2 = Parameter(name='outer acyl chain thickness', value=12.0).range(8, 16)
vf_overlayer = Parameter(name='volume fraction overlayer', value=0.05).range(0.0, 0.2)
l_lipid_overlayer = Parameter(name='overlayer acyl chain thickness', value=10.0).range(8, 25)
dz_overlayer = Parameter(name='overlayer acyl chain separation distance', value=20.0).range(15.0, 50.0)
sigma = Parameter(name='bilayer roughness', value=3).range(0.5, 9)
global_rough = Parameter(name ='substrate roughness', value=5).range(2, 9)
d_oxide = Parameter(name='silicon oxide layer thickness', value=10).range(5, 30)
d_Cr =  Parameter(name='chromium layer thickness', value=40).range(10, 150)
d_gold =  Parameter(name='gold layer thickness', value=140).range(100, 200) #thickness of gold
rough_cr_au =  Parameter(name='gold chromium roughness', value=10).range(2, 24.0) # roughness of Cr/Au interface
nf_tether =  Parameter(name='number fraction tether', value=0.45).range(0.2, 1.0) # number fraction of tether molecules in inner leaflet
mult_tether =  Parameter(name='bME to tether ratio', value=3).range(0.1, 4) #ratio of bME to tether molecules at surface
l_tether =  Parameter(name='tether length', value=10).range(3, 18) #distance from substrate to inner headgroup/acyl chain interface

### Define bilayer object
DOPC = cmp.Lipid(name='DOPC', headgroup=cmp.pc, tails=2 * [cmp.oleoyl], methyls=[cmp.methyl])
blm = mol.tBLM(tether=lipids.HC18SAc, filler=cmp.bme, lipids=[DOPC, lipids.DOPE],
                lipid_nf=[0.5, 0.5])

betap = Parameter(name='beta rotation', value=110.0).range(0, 180.)
gammap = Parameter(name='gamma rotation', value=10.0).range(0, 360.)
nf_prot = Parameter(name='number fraction protein', value=0.003).range(0.0, 0.02)
protz = Parameter(name='protein center of mass position', value=100.0).range(50, 150)
protexchratio = Parameter(name='proton exchange efficiency', value=0.9)

vf_bilayer_prot = Parameter(name='volume fraction bilayer protein', value=0.95).range(0.0, 1.0)
#vf_bilayer_prot = vf_bilayer
dl_lipid2 = Parameter(name='outer acyl chain thickness change', value=0.0).range(-3, 3)
radius_defects = 10.**Parameter(name='log10 radius_defects', value=8.0) #.range(0, 8)

l_mr = Parameter(name='length missing residues', value=10).range(8.0, 50.0)
frac_pos_mr = Parameter(name='fractional position missing residues', value=0.5).range(0.0, 1.0)

# define protein object: Tubulin PDB with missing residues
# prot = mol.ContinuousEuler(mol.pdbto8col('1TVK_no_hetatm.pdb', 'tubulin8col.txt'), name='tubulin_1TVK')
prot = mol.ContinuousEuler('tubulin8col.txt', name='tubulin_1TVK')
protm = mol.ContinuousEulerMissingResidues(prot, name='Tubulin')
mr = protm.add_missing_residues('QMPSDKTIGGGDDSFNTFFSETGAGK', (34, 61))
mr.name='tubulin_missing_residues'
blmprot = mol.BLMProteinComplex(blms=[blm], proteins=[protm])

## === Stack ===
##
## First, we create a 'material' for each bulk layer, which has an real and imaginary
## scattering length density, stored in a Refl1d object called 'SLD'
d2o = SLD(name='D2O buffer', rho=6.3000, irho=0.0000)
h2o = SLD(name='H2O buffer', rho=-0.56, irho=0.0000)
d2o_prot = SLD(name='D2O buffer with protein', rho=6.3000, irho=0.0000)
h2o_prot = SLD(name='H2O buffer with protein', rho=-0.56, irho=0.0000)
air = SLD(name='air', rho=0.000, irho=0.0000)
siox = SLD(name='silicon oxide', rho=4.1000, irho=0.0000)
silicon = SLD(name='silicon', rho=2.0690, irho=0.0000)
cr = SLD(name='chromium', rho=2.7, irho=0.0)
gold = SLD(name='gold', rho=4.4, irho=0.0) #iro is the absorption of neutrons, should be 0

## Then bulk layers are created, each with its own 'material'.  If you want to force
## two layers to always match SLD you can use the same material in multiple layers.
## The roughnesses of each layer are set to zero to begin with:

layer_d2o = Slab(material=d2o, thickness=0.0000, interface=5.0000)
layer_h2o = Slab(material=h2o, thickness=0.0000, interface=5.0000)
layer_d2o_prot = Slab(material=d2o_prot, thickness=0.0000, interface=5.0000)
layer_h2o_prot = Slab(material=h2o_prot, thickness=0.0000, interface=5.0000)
layer_air = Slab(material=air, thickness=0.0000, interface=5.0000)
layer_siox = Slab(material=siox, thickness=d_oxide, interface=global_rough)
layer_silicon = Slab(material=silicon, thickness=0.0000, interface=global_rough)
layer_cr = Slab(material=cr, thickness=d_Cr, interface=rough_cr_au)
layer_gold = Slab(material=gold, thickness=d_gold - (blm.substrate.z + 0.5 * blm.substrate.length), interface=0.0000)

substrate = layer_silicon  | layer_siox | layer_cr | layer_gold

## Use the bilayer definition function to generate the bilayer SLD profile, passing in the relevant parameters.
## Note that substrate and bulk SLDs are linked to their respective materials.

sample, sampleh = make_samples(bilayer, substrate, [d2o, h2o], sigma=sigma, 
                             global_rough=global_rough, rho_substrate=gold.rho,
                             l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2,
                             nf_tether=nf_tether, mult_tether=mult_tether,
                             vf_bilayer=vf_bilayer)

sample_prot, sample_proth = make_samples(bilayer_prot, substrate, [d2o_prot, h2o_prot], sigma=sigma, 
                             global_rough=global_rough, rho_substrate=gold.rho,
                             l_tether=l_tether, l_lipid1=l_lipid1, l_lipid2=l_lipid2 + dl_lipid2,
                             nf_tether=nf_tether, mult_tether=mult_tether,
                             vf_bilayer=vf_bilayer_prot, gamma=gammap, beta=betap, protz=protz, protnf=nf_prot,
                             radius_defects=radius_defects, l_mr=l_mr, frac_pos_mr=frac_pos_mr)

## Set sample parameter ranges and constraints between layer properties, if these are not set using parameters previously

# nSLD parameters
d2o.rho.range(5.3000, 6.5000)
d2o_prot.rho.range(5.3000, 6.5000)
h2o.rho.range(-0.6, 1.6)
h2o_prot.rho.range(-0.6, 1.6)
siox.rho.range(3.6000, 5.1000)
cr.rho.range(2.7000, 3.6000)
gold.rho.range(4.2000, 4.8000)

## === Data files ===

probe = load4('mb197_d2o.ort', back_reflectivity=True, name='D2O buffer')
probeh = load4('mb198_h2o.ort', back_reflectivity=True, name='H2O buffer')
probe_prot = load4('mb203_d2o_tubulin.ort', back_reflectivity=True, name='D2O + tubulin')
probe_proth = load4('mb205_h2o_tubulin.ort', back_reflectivity=True, name='H2O + tubulin')

# Background parameter
# probe.background.value = 0.0000
probe.background.range(-1e-7, 1e-5)
probeh.background.range(-1e-7, 1e-5)
probe_prot.background.range(-1e-7, 1e-5)
probe_proth.background.range(-1e-7, 1e-5)
probe.intensity.range(0.90, 1.15)
probe_proth.intensity = probe_prot.intensity = probeh.intensity = probe.intensity
probe.theta_offset.range(-0.02, 0.02)
probe_prot.theta_offset = probe_proth.theta_offset = probeh.theta_offset = probe.theta_offset
probe.sample_broadening.range(-0.003, 0.02)
probe_prot.sample_broadening = probe_proth.sample_broadening = probeh.sample_broadening = probe.sample_broadening

# Define critical edge oversampling for samples that require it (typically D2O only)
probe.critical_edge(substrate=silicon, surface=d2o)
probe_prot.critical_edge(substrate=silicon, surface=d2o_prot)

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
model_prot = Experiment(sample=sample_prot, probe=probe_prot, dz=STEPSIZE, step_interfaces = step)
model_proth = Experiment(sample=sample_proth, probe=probe_proth, dz=STEPSIZE, step_interfaces = step)
problem = FitProblem([model, modelh, model_prot, model_proth])

problem.blms = [blm]
problem.proteins = [protm]
problem.blmprot = blmprot

def custom_plot():

    import plotly.graph_objs as go
    from refl1d.webview.server.colors import COLORS

    moldat = problem.moldat

    def hex_to_rgb(hex_string):
        r_hex = hex_string[1:3]
        g_hex = hex_string[3:5]
        b_hex = hex_string[5:7]
        return int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)

    n_lipids = 2
    group_names = {'gold substrate': ['bilayer.substrate'],
               #'silicon oxide': ['bilayer.siox'],
               'bME': ['bilayer.bME'],
               'tether': ['bilayer.tether_bme', 'bilayer.tether_free', 'bilayer.tether_hg'],
               'tether acyl chains': ['bilayer.tether_methylene', 'bilayer.tether_methyl'],
               'inner headgroups': [f'bilayer.headgroup1_{i}' for i in range(1, n_lipids + 1)],
               'inner acyl chains': [f'bilayer.methylene1_{i}' for i in range(1, n_lipids + 1)] + [f'bilayer.methyl1_{i}' for i in range(1, n_lipids + 1)],
               'outer acyl chains': [f'bilayer.methylene2_{i}' for i in range(1, n_lipids + 1)] + [f'bilayer.methyl2_{i}' for i in range(1, n_lipids + 1)],
               'outer headgroups': [f'bilayer.headgroup2_{i}' for i in range(1, n_lipids + 1)]
              }
    
    overlay_groups = {'tubulin': ['protein.tubulin_1TVK'],
                      'missing residues': ['protein.tubulin_missing_residues']}
    #overlay_groups = {}
    
    normarea = moldat['bilayer.normarea']['area']

    fig = go.Figure()
    traces = []
    MOD_COLORS = COLORS[1:]
    color_idx = 1
    sumarea = 0
    sumarea_overlay = 0
    area_gps = {}
    for lbl, item in group_names.items():
        area = 0
        for gp in item:
            if gp in moldat.keys():
                zaxis = moldat[gp]['zaxis']
                area += moldat[gp]['area']
            else:
                print(f'Warning: {gp} not found')
        area_gps[lbl] = area
        sumarea += area

    overlay_gps = {}
    for lbl, item in overlay_groups.items():
        area = 0
        for gp in item:
            if gp in moldat.keys():
                zaxis = moldat[gp]['zaxis']
                area += moldat[gp]['area']
            else:
                print(f'Warning: {gp} not found')
        overlay_gps[lbl] = area
        sumarea_overlay += area

    frac_replacement = np.ones_like(normarea)
    frac_replacement[(sumarea + sumarea_overlay) > normarea] = (sumarea / (normarea - sumarea_overlay))[(sumarea + sumarea_overlay) > normarea]

    for lbl, area in area_gps.items():
        color = MOD_COLORS[color_idx % len(MOD_COLORS)]
        plotly_color = ','.join(map(str, hex_to_rgb(color)))
        traces.append(go.Scatter(x=zaxis,
                                    y=area / normarea / frac_replacement,
                                    mode='lines',
                                    name=lbl,
                                    line=dict(color=color)))
        traces.append(go.Scatter(x=zaxis,
                                    y=area / normarea / frac_replacement,
                                    mode='lines',
                                    line=dict(width=0),
                                    fill='tozeroy',
                                    fillcolor=f'rgba({plotly_color},0.3)',
                                    showlegend=False
                                    ))
        color_idx += 1

    for lbl, area in overlay_gps.items():
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
    problem.custom_plot().show()

