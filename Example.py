#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  21 11:13:35 2020

@author: sam
PDE/Wave Equation solver
"""
import numpy as np
import pylab as plt

from WaveEquation import WaveSolver

# Creating arrays for domain and displacement
dx = 20e-6
dy = 20e-6
x = np.arange(0, 20e-3, dx)
y = np.arange(0, 20e-3, dy)

# Wavespeed
c = np.ones([len(x), len(y)])*3000
c[500:, 500:] = 1500

# Applying phase delay to excitation signal
angle = 25 * np.pi/180
Y, X = np.meshgrid(y, x)
R = np.sqrt((X - X[750, 500])**2 + (Y - Y[750, 500])**2)
t0 = -(R - R[500, 1])/3000. + 5e-6

A = np.zeros(c.shape)
A[100:900, 1] = 1
# Create wavesolver instance
Soln = WaveSolver(x, y, c, dt=2.5e-9, params=(A, 2e6, 4e-7, t0))

# Get first timestep to set up plot
u = Soln.get_snapshot()

# Create image plot
Extent = (1e3*x.min(), 1e3*x.max(), 1e3*y.min(), 1e3*y.max())
fig, ax = plt.subplots()
im = ax.imshow(u.T, aspect='auto', origin='lower', extent=Extent)
im.set_clim(-4.5, 4.5)
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')
plt.tight_layout()

for n in range(8000):
    # Iterate forward by 20 us (8000*2.5ns), showing every 100th timestep
    Soln.solve_step()
    if not n%100:
        u = Soln.get_snapshot()
        im.set_data(u.T)
        plt.pause(0.01)
