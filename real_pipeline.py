"""
BAO Project - Submission date 14th September 2026

This is the code used to run the real dataset. 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from nbodykit.lab import *
from nbodykit import setup_logging
from astropy.io import fits
from nbodykit.lab import ArrayCatalog

setup_logging()
# 1. Load the real dataset also using a subsample of it

hdul = fits.open('/home/pasta/bao_project/eBOSS_LRG_clustering_data-NGC-vDR16.fits')
hdul_ran = fits.open('/home/pasta/bao_project/eBOSS_LRG_clustering_random-NGC-vDR16.fits')

idx = np.random.choice(5460719, 2150000, replace=False)

data = hdul[1].data
data_ran = hdul_ran[1].data

data_ran_sub = data_ran[idx]

weight_total = data['WEIGHT_SYSTOT'] * data['WEIGHT_CP'] * data['WEIGHT_NOZ'] * data['WEIGHT_FKP']
weight_ran = data_ran_sub['WEIGHT_SYSTOT'] * data_ran_sub['WEIGHT_CP'] * data_ran_sub['WEIGHT_NOZ'] * data_ran_sub['WEIGHT_FKP']

# 2. Cosmology + linear power spectrum + comoving distance calculation

cosmos = cosmology.Planck15
Plin = cosmology.LinearPower(cosmos, redshift= 0.55, transfer='EisensteinHu')
Plin_nowiggle = cosmology.LinearPower(cosmos, redshift = 0.55, transfer = 'NoWiggleEisensteinHu')
h = cosmos.h

cf_nowiggle = cosmology.CorrelationFunction(Plin_nowiggle)

# Values of comving distance is already in Mpc. 
dis = cosmos.comoving_distance(data['Z'])
dis_ran = cosmos.comoving_distance(data_ran_sub['Z'])

# 3. Wrapping catalogs and attaching weights

data_cat = ArrayCatalog(data)
random_cat = ArrayCatalog(data_ran_sub)
data_cat['WEIGHT_TOTAL'] = weight_total
random_cat['WEIGHT_TOTAL'] = weight_ran

# 4. Two-point correlation function (periodic box case)

r_edge = np.linspace(1, 200, 20)
twopcf = SurveyData2PCF(mode='1d', data1=data_cat, randoms1=random_cat, edges=r_edge, cosmo=cosmos, ra='RA', dec='DEC', redshift='Z', weight='WEIGHT_TOTAL')

r = twopcf.corr['r']
xi = twopcf.corr['corr']

xi_scaled = xi * r**2

np.save('r_values.npy', r)
np.save('xi_values.npy', xi) #just to save on computing time.

#Plotting 

r_hunits = r * h
xi_nowiggle = cf_nowiggle(r_hunits)

plt.plot(r, r**2 * xi, label='Real Data')
plt.plot(r, r**2 * xi_nowiggle, label='No-Wiggle Theory', linestyle='--')
plt.axvspan(80, 110, alpha=0.3, color='blue', label='possible feature (Region 1)')
plt.axvspan(145, 160, alpha=0.3, color='green', label='Possible feature (Region 2)')
plt.axvline(150, color='black', linestyle=':', label='Expected BAO scale (~150 Mpc)')
plt.xlabel(r'$r$ [Mpc]')
plt.ylabel(r'$r^2 \xi(r)$')
plt.title('Real Data vs. No-Wiggle Theory: Searching for BAO Bump')
plt.legend()
plt.savefig('bao_comparison_plot.png')
plt.show()