#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 09:22:09 2013

@author: pete
"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np 
import foampy
from gendynmeshdict import meantsr

t, torque, drag = foampy.load_all_torque_drag(torque_axis="x")

# Compute tip speed ratio
R = 0.13
U = 1.0
omega = -meantsr*U/R 		# direction is reversed
theta = omega*t*180.0/np.pi

print("Mean tsr:", meantsr)

# Pick an index to start from for mean calculations and plotting
# (allow turbine to reach steady state)
try:
    i = np.where(np.round(theta) == 360)[0][0]
except IndexError:
    i = 5
i2 = -1

# Compute power coefficient
area = np.pi*R**2
power = torque*omega
cp = power/(0.5*1000*area*U**3)
print("Mean cp:", np.mean(cp[i:i2]))

# Compute drag coefficient
cd = drag/(0.5*1000*area*1**2)
print("Mean cd:", np.mean(cd[i:i2]))

plt.close('all')
plt.plot(theta[i:i2], cp[i:i2])
plt.title(r"$\lambda = %1.1f$" %meantsr)
plt.xlabel(r"$\theta$ (degrees)")
plt.ylabel(r"$C_P$")
#plt.ylim((0, 0.7))
plt.show()
