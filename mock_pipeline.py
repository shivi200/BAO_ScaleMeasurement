"""
BAO Project - Submission date 14th September 2026

Purpose: build and validate the nbodykit correlation-function + power-spectrum
pipeline on a synthetic lognormal mock catalog with known BAO wiggles,
before the real SDSS/DESI catalog.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from nbodykit.lab import *
from nbodykit import setup_logging

setup_logging()

# 1. Cosmology + linear power spectrum

cosmo = cosmology.Planck15
Plin = cosmology.LinearPower(cosmo, redshift= 0.55, transfer='EisensteinHu')
""" Computes the linear power spectrum using Planck15 dataset, and Eisenstein
&Hu transfer function that doesn't use Boltzmann solvers (slow). 
"""

# 2. Generate a lognormal mock catalog with BAO wiggles baked in

BoxSize = 1380.0
nbar = 3e-4
bias = 2.0

cat = LogNormalCatalog(Plin=Plin, nbar=nbar, BoxSize=BoxSize,Nmesh=256, bias=bias, seed=42)

print(f"Mock catalog created: {cat.csize} objects in a "f"({BoxSize} Mpc/h)^3 periodic box")

# 3. Two-point correlation function (periodic box case)

r_edges = np.linspace(1, 200, 50)

twopcf = SimulationBox2PCF(mode='1d', data1=cat, edges=r_edges, BoxSize=BoxSize, periodic=True, nthreads=10)

r = twopcf.corr['r']
xi = twopcf.corr['corr']


# 4. Power spectrum via FFT

mesh = cat.to_mesh(Nmesh=256, BoxSize=BoxSize, resampler='tsc',compensated=True)
r_power = FFTPower(mesh, mode='1d', dk=0.005, kmin=0.01)

k = r_power.power['k']
Pk = r_power.power['power'].real
shotnoise = r_power.power.attrs.get('shotnoise', 0.0)
Pk_shot_sub = Pk - shotnoise

# 5. Generate a no-wiggle graph for comparison 

Plin_nowiggle = cosmology.LinearPower(cosmo, redshift = 0.55, transfer = "NoWiggleEisensteinHu")
Pk_nowiggle = Plin_nowiggle(k) * bias**2
ratio = Pk_shot_sub / Pk_nowiggle
broadband = savgol_filter(ratio, window_length = 41, polyorder = 3)
wiggle = ratio / broadband


# 6. Plot both diagnostics

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].loglog(k, Pk_shot_sub, label = "Mock (with Wiggle)")
axes[0].loglog(k, Pk_nowiggle, label = "No-Wiggle reference", linestyle = '--')
axes[0].set_xlabel(r'$k$ [$h$/Mpc]')
axes[0].set_ylabel(r'$P(k) [(Mpc/$h)^3$]')
axes[0].legend()
axes[0].set_title('Power spectrum: mock vs no-wiggle')

axes[1].plot(r, r**2 * xi)
axes[1].axvline(105, color='gray', ls='--', label='approx. BAO scale')
axes[1].set_xlabel(r'$r$  [Mpc/$h$]')
axes[1].set_ylabel(r'$r^2\,\xi(r)$')
axes[1].legend()
axes[1].set_title('Two-point correlation function (mock)')

axes[2].plot(k, ratio, label = "Ratio without filter", linestyle = 'dotted')
axes[2].plot(k, wiggle, label = "Ratio with filter")
axes[2].axhline(1.0, color = 'red', ls = '--')
axes[2].set_xlabel(r'$k$ [$h$/Mpc]')
axes[2].set_ylabel(r'$Ratio$')
axes[2].legend()
axes[2].set_title('Mock / no wiggle (isloated BAO wiggle)')

plt.tight_layout()
plt.savefig('bao_mock_diagnostics.png', dpi=150)
plt.show()

print("Pipeline complete. Check bao_mock_diagnostics.png:")