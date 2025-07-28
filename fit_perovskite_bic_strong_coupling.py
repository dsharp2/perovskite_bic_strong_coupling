import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

e_exc = 2.283135726 # Define exciton energy in eV. Example values are given here
gam_cav = 0.001254334 # Define caviy half-width half maximum (HWHM) in eV
gam_exc = 0.0172 # Define exciton HWHM in eV

m_0 = 1.3205111040631055 # Define linear guided mode dispersion from simulation
b_0 = 2.497169985164217

# Define lower polariton energy and return real part
def polariton(k, g, m, b):
	return np.real((1/2) * (e_exc + m*k + b + (gam_cav + gam_exc)*1j) - np.sqrt(g**2 + (1/4) * (e_exc - (m*k + b) + (gam_exc - gam_cav)*1j)**2))

# Define guided mode dispersion
def gmr(k, m, b):
	return m*k+b

# Import wavevectors (k) and lower polariton energies
k_data = np.genfromtxt('...csv', skip_header=0, skip_footer=0, dtype=None,usecols=0,delimiter=',',encoding=None).astype(float)

lp_data = np.genfromtxt('....csv', skip_header=0, skip_footer=0, dtype=None,usecols=1,delimiter=',',encoding=None).astype(float)

# Fit experimental lower polariton data
popt, pcov = curve_fit(polariton, k_data, lp_data, bounds = ([1e-3, m_0*0.9, b_0*0.9], [50e-3, m_0*1.1, b_0*1.1]))

# Return fitted parameters and plot
[g_fit, m_fit, b_fit] = popt

plt.plot(k_data, lp_data, c='k', label = 'LP Data')
plt.scatter(k_data, polariton(k_data, g_fit, m_fit, b_fit), color='k', label = 'LP Fit')
plt.plot(k_data, gmr(k_data, m_fit, b_fit), c='b', linestyle='--', label = 'GMR')
plt.axhline(e_exc, c='r', label='Exciton')
plt.legend()