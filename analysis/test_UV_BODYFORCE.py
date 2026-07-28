import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from helpers import prep_ds, load_yaml

def analytical_u(F, t):
    return F*t


experimentpath = "runs/test_UV_BODYFORCE"
params = load_yaml(f"{experimentpath}/resolved_config.yaml")

# Analytical solution for u
NTIMES = params["time_stepping"]["NTIMES"]
DT = params["time_stepping"]["DT"]
NHIS = params["time_stepping"]["NHIS"]
F = params["bodyforce"]["BFRC_U"]*1e-7 

dt = DT*NHIS
T = NTIMES*DT
t = np.arange(0, T+dt/2, dt)
u_analytical = analytical_u(F, t)

# Read the dataset and prepare it
time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
ds = xr.open_dataset(params["io"]["output_dir"]+"/"+params["files"]["his"], decode_times=time_coder)
ds, grid = prep_ds(ds, params)
u_simulation = ds.u.values


# check that the analytical solution matches the numerical solution for every xi, eta and s_rho point
# Print message if the assertion is successful
for i in range(ds.xi_u.size): 
    for j in range(ds.eta_rho.size):
        for k in range(ds.s_rho.size):
            assert np.isclose(u_analytical, u_simulation[:, k, j, i]).all(), f"Mismatch at xi={i}, eta={j}, s_rho={k}"
print("Simulated u is equal to analytical u for all grid points.")


# Also plot the analytical and numerical solutions for visual inspection
fig, ax = plt.subplots(figsize=(10, 6))
th = t/(60*60)

# Plot the analytical solution
ax.plot(th, u_analytical, label='Analytical Solution', color='blue', linewidth=2)

# Plot the numerical solution with spread
u_simulation_mean = np.mean(u_simulation, axis=(1, 2, 3))
ax.plot(th, u_simulation_mean, 
        label='Numerical solution (averaged over all points)', color='orange', linestyle='--', linewidth=2)
u_simulation_std = np.std(u_simulation, axis=(1, 2, 3))
ax.fill_between(th, u_simulation_mean - u_simulation_std, u_simulation_mean + u_simulation_std, 
                color='orange', alpha=0.2, label='Numerical solution spread (1 std dev)')

ax.set_xlabel('Time (hours)')
ax.set_ylabel('u (m/s)')
ax.set_title('Comparison of analytical and numerical solutions for u given a constant body force')
ax.legend()
ax.grid()

# add anlytical expression for u in the plot
ax.text(0.5, 0.1, r'$u(t) = F t$',
        transform=ax.transAxes, fontsize=12, verticalalignment='bottom', 
        horizontalalignment='center')

fig.savefig(f"figures/test_UV_BODYFORCE.png")