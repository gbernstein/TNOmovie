#!/usr/bin/env python
# coding: utf-8

# In[3]:


import matplotlib.pyplot as pl
from astropy.io import fits
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpatches
from astropy import wcs
from astropy.wcs import WCS
import astropy
from astropy.utils.data import get_pkg_data_filename
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, get_body_barycentric_posvel, get_body, get_body_barycentric, SkyCoord
from astropy import units as u
import numpy as np
import scipy
import pandas as pd

SunGM       = 4.*np.pi*np.pi/1.000037773533  #solar gravitation
MercuryGM   = 6.55371264e-06
VenusGM     = 9.66331433e-05
EarthMoonGM = 1.20026937e-04
MarsGM      = 1.27397978e-05;
JupiterGM   = 3.76844407e-02 + 7.80e-6
SaturnGM    = 1.12830982e-02 + 2.79e-6
UranusGM    = 1.72348553e-03 + 0.18e-6
NeptuneGM   = 2.03318556e-03 + 0.43e-6
SolarSystemGM = SunGM + MercuryGM + VenusGM + EarthMoonGM + MarsGM + JupiterGM + SaturnGM + UranusGM + NeptuneGM
EclipticInclination = 23.43928 * np.pi/180

import numba
@numba.vectorize(["float32(float32,float32)", "float64(float64,float64)"],nopython=True)
def solve_anomaly(e, M0):
	#Args:
		#e (float, array): eccentricity
		#M0 (float, array): Mean anomaly
	#Returns:
		#E (float, array): Eccentric anomaly at input epoch (i.e. in M0)
    sol = M0
    delta = 0.0
    ones = 1.0
    for i in range(1000):
        delta = (M0 - (sol - e * np.sin(sol)))/(ones - e*np.cos(sol))
        sol += delta
    return sol

def q(E, a, e):
    q1 = a*(np.cos(E)-e)
    q2 = a*np.sqrt(1-e**2)*np.sin(E)
    if np.isscalar(a):
        return np.array([q1,q2,0])
    else:
        return np.array([q1,q2,np.zeros_like(a)])

def R(omega, Omega, i):
    cO = np.cos(np.pi*Omega/180)
    sO = np.sin(np.pi*Omega/180)
    co = np.cos(np.pi*omega/180)
    so = np.sin(np.pi*omega/180)
    ci = np.cos(np.pi*i/180)
    si = np.sin(np.pi*i/180)
    R = np.array([[cO * co - sO * so * ci, - cO * so - sO * co * ci, sO * si],
                [sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
                [si * so, si*co, ci]])
    if np.isscalar(omega):
        return R
    else:
        return np.transpose(R,(2,0,1))

def dqdt(E, a, e, mu):
    n = np.sqrt(mu)/np.power(a,3./2)
    den = 1 - e*np.cos(E)
    q1 = -n*a*np.sin(E)/den
    q2 = n*a*np.sqrt(1-e**2)*np.cos(E)/den
    if np.isscalar(a):
        return np.array([q1,q2,0])
    else:
        return np.array([q1,q2,np.zeros_like(a)])
    
def rotate_to_ecliptic(xv, inverse = False):
    cosEcl = np.cos(EclipticInclination)
    sinEcl = -np.sin(EclipticInclination)
    if inverse:
        sinEcl = - sinEcl
    xv[:,1], xv[:,2] = cosEcl * xv[:,1] + sinEcl * xv[:,2], -sinEcl * xv[:,1] + cosEcl * xv[:,2]
    xv[:,4], xv[:,5] = cosEcl * xv[:,4] + sinEcl * xv[:,5], -sinEcl * xv[:,4] + cosEcl * xv[:,5]
    return xv

def keplerian_to_cartesian(keplerian, epoch, helio = False, ecliptic = False):
    mu = SunGM if helio else SolarSystemGM
    M0 = np.array((epoch - keplerian[:,5])*(np.sqrt(mu)/np.power(keplerian[:,0],3./2)))
    E = solve_anomaly(np.array(keplerian[:,1]), M0)
    q_vec = q(E, keplerian[:,0], keplerian[:,1])
    Rot = R(keplerian[:,4], keplerian[:,3], keplerian[:,2])
    x = np.einsum('...ij,j...', Rot, q_vec)
    dqdt_vec = dqdt(E, keplerian[:,0], keplerian[:,1], mu)
    v = np.einsum('...ij,j...', Rot, dqdt_vec)
    xv = np.zeros_like(keplerian)
    xv[:,0:3] = x
    xv[:,3:] = v
    if not ecliptic:
        xv = rotate_to_ecliptic(xv)
    return xv

def cartesian_to_keplerian(cartesian, epoch, helio = False, ecliptic = False):
    mu = SunGM if helio else SolarSystemGM
    xv = np.zeros_like(cartesian)
    if not ecliptic:
        cosEcl = np.cos(EclipticInclination)
        sinEcl = np.sin(EclipticInclination)
        xv[:,0] = cartesian[:,0]
        xv[:,3] = cartesian[:,3]
        xv[:,1], xv[:,2] = cosEcl * cartesian[:,1] + sinEcl * cartesian[:,2], -sinEcl * cartesian[:,1] + cosEcl * cartesian[:,2]
        xv[:,4], xv[:,5] = cosEcl * cartesian[:,4] + sinEcl * cartesian[:,5], -sinEcl * cartesian[:,4] + cosEcl * cartesian[:,5]
    else:
        xv = cartesian
    x = np.sqrt(xv[:,0]*xv[:,0] + xv[:,1]*xv[:,1]+ xv[:,2]*xv[:,2])
    vsq_mu = (xv[:,3]**2 + xv[:,4]**2 + xv[:,5]**2)/mu
    inv_a = 2./x - vsq_mu
    a = 1./inv_a
    x_dot_v = (xv[:,0] * xv[:,3] + xv[:,1] * xv[:,4] + xv[:,2] * xv[:,5])
    pref = (vsq_mu - 1./x)
    e_vec = (pref * xv[:,0:3].T -  x_dot_v * xv[:,3:6].T/mu).T
    e = np.sqrt(e_vec[:,0] * e_vec[:,0] + e_vec[:,1] * e_vec[:,1] + e_vec[:,2] * e_vec[:,2])
    h_vec = np.cross(xv[:,0:3], xv[:,3:6])
    n_vec = np.cross(np.array([0,0,1]), h_vec)
    h = np.sqrt(h_vec[:,0] * h_vec[:,0] + h_vec[:,1] * h_vec[:,1] + h_vec[:,2] * h_vec[:,2])
    n = np.sqrt(n_vec[:,0] * n_vec[:,0] + n_vec[:,1] * n_vec[:,1] + n_vec[:,2] * n_vec[:,2])
    cos_i = h_vec[:,2]/h 
    cosOmega = n_vec[:,0]/n 
    cosomega = (n_vec[:,0] * e_vec[:,0] + n_vec[:,1] * e_vec[:,1] + n_vec[:,2] * e_vec[:,2])/(n*e)
    i = np.arccos(cos_i) * 180./np.pi
    Omega = np.arccos(cosOmega) * 180./np.pi
    omega = np.arccos(cosomega) * 180./np.pi
    Omega[np.where(n_vec[:,1] < 0)] = 360 - Omega[np.where(n_vec[:,1] < 0)]
    omega[np.where(e_vec[:,2] < 0)] = 360 - omega[np.where(e_vec[:,2] < 0)]
    p = h*h/mu
    b = a * np.sqrt(1 - e*e)
    xbar = (p - x)/e 
    ybar = x_dot_v/e * np.sqrt(p/mu)
    E = np.arctan2(ybar/b, xbar/a + e)
    M = E - e*np.sin(E)
    T_p = epoch - M * np.sqrt(a**3/mu)
    aei = np.zeros_like(xv)
    aei[:,0] = a
    aei[:,1] = e 
    aei[:,2] = i 
    aei[:,3] = Omega
    aei[:,4] = omega
    aei[:,5] = T_p
    return aei


# In[ ]:




