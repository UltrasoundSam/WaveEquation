#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  21 11:13:35 2020

@author: sam
PDE/Wave Equation solver
"""
import numpy as np
import pylab as plt

from waveequation import WaveSolver


def main():
    # Creating arrays for domain and displacement
    dx = 20e-6
    dy = 20e-6
    x = np.arange(0, 20e-3, dx)
    y = np.arange(0, 20e-3, dy)

    # Wavespeed
    c = np.ones([len(x), len(y)])*3000
    c[500:, 500:] = 1500

    # Applying phase delay to excitation signal
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
    im.set_clim(-7e-17, 7e-17)
    ax.set_xlabel('X (mm)')
    ax.set_ylabel('Y (mm)')
    ax.set_title(r'Time - 0.00 $\mu$s')
    plt.tight_layout()

    max_val = 0
    for n in range(8001):
        # Iterate forward by 20 us (8000*2.5ns), showing every 100th timestep
        Soln.solve_step()
        if not n % 100:
            u = Soln.get_snapshot()
            im.set_data(u.T)
            ax.set_title(r'Time - {0:0.2f} $\mu$s'.format(2.5e-3*n))
            if (newmax := np.abs(u.max())) > max_val:
                max_val = newmax
            plt.pause(0.01)


if __name__ == '__main__':
    main()
