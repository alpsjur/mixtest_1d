'''
The idea if this test is to see if the production term
in STRUCTURE_MIXING balance dissipation in steady state,
when no other production term is present (i.e. no shear, 
no stratification).    
'''

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from helpers import prep_ds, load_yaml

def calc_Pd(str_a, Cd, u):
    return 0.5*Cd*str_a*np.abs(u**3)


experimentpath = "runs/test_STRUCTURE_PRODUCTION"
params = load_yaml(f"{experimentpath}/resolved_config.yaml")

str_a = params["structure"]["str_a"]
Cd = params["structure"]["CD"]

# Read the dataset and prepare it
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
ds = xr.open_dataset(params["io"]["output_dir"]+"/"+params["files"]["his"], decode_times=time_coder)
ds, grid = prep_ds(ds, params)


# shift all the variables to the cell centers for easier comparison
u = grid.interp(ds.u, 'X')
tke = grid.interp(ds.tke, 'Z')
gls = grid.interp(ds.gls, 'Z')
Pd = calc_Pd(str_a, Cd, u)

TKE = grid.integrate(tke, ('X', 'Y', 'Z'))
GLS = grid.integrate(gls, ('X', 'Y', 'Z'))
PD = grid.integrate(Pd, ('X', 'Y', 'Z'))    

fig, ax = plt.subplots(figsize=(6, 4))

ax.plot(PD, label=r'Analytical $P_d$', color='blue')
ax.plot(GLS, label=r'Numerical $\epsilon$', color='orange', linestyle='--')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Domain integrated values')
ax.legend()
ax.set_title('Comparing domain integrated structure production\nand dissipation terms for shear free, unstratified flow')
ax.grid()

fig.savefig(f"figures/test_STRUCTURE_PRODUCTION.png")