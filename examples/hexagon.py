#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thurs Nov  18 11:13:35 2021

@author: sam
PDE/Wave Equation solver
"""
import numpy as np
import pylab as plt

from matplotlib.patches import RegularPolygon
from matplotlib.animation import FuncAnimation

from waveequation import WaveSolver


def init(fig, ax, im):
    '''
    Initiate the limits and labels of the plot
    '''
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_xlabel(r'X (mm)')
    ax.set_ylabel(r'Y (mm)')
    ax.set_title(r'Time - 0.00 $\mu$s')
    im.set_clim(-1, 1)
    cbar = fig.colorbar(im)
    cbar.set_ticks([])
    fig.tight_layout()
    return im,


def animate(i, im, ax, WaveSolve):
    '''
    Advance the wavesolver by 0.1 microseconds and then update plot.
    '''
    # Update plot by 0.05 microseconds (20 iterations)
    for _ in range(20):
        WaveSolve.solve_step()
    u = WaveSolve.get_snapshot()

    # Update plot data and titles
    im.set_data(u.T)
    ax.set_title(r'Time - {0:0.2f} $\mu$s'.format(1e6*WaveSolve.t))
    max_val = np.abs(u).max()
    im.set_clim(-max_val, max_val)
    return im,


def hexagon_mask(x, y, radius):
    '''
    Creates hexagon-shaped mask for defining the wave speed in a
    hexagonal shape. Defines hexagon with a given radius.
    '''
    # Create meshgrid
    X, Y = np.meshgrid(np.abs(x), np.abs(y))

    # Perform two checks to see if you are inside hexagon
    check1 = ((X <= radius)) & (Y <= radius*np.sqrt(3)/2)
    check2 = (Y <= np.sqrt(3) * (radius - X))

    # Combine two checks to create mask
    checkt = check1 & check2
    return checkt


def main():
    # Creating arrays for domain and displacement
    dx = 20e-6
    dy = 20e-6
    x = np.arange(-10e-3, 10e-3, dx)
    y = np.arange(-10e-3, 10e-3, dy)

    # Wavespeed - define a hexagon
    c = np.ones([len(x), len(y)])*10
    radius = 10e-3

    mask = hexagon_mask(x, y, radius)
    c[mask] = 3000

    A = np.zeros(c.shape)
    A[499:501, 499:501] = 1

    # Add hexagon on plot to delineate the domain
    hexagon = RegularPolygon((0, 0), numVertices=6,
                             radius=1e3*radius, alpha=.5)

    # Create wavesolver instance
    Soln = WaveSolver(x, y, c, dt=2.5e-9, params=(A, 2e6, 4e-7, 0e-3))

    # Get first timestep to set up plot
    u = Soln.get_snapshot()

    # Create image plot
    Extent = (1e3*x.min(), 1e3*x.max(), 1e3*y.min(), 1e3*y.max())
    fig, ax = plt.subplots()
    ax.add_patch(hexagon)
    im = ax.imshow(u.T, aspect='auto', origin='lower',
                   extent=Extent, clip_path=hexagon, clip_on=True)

    init(fig, ax, im)

    # Animate through
    anim = FuncAnimation(fig, animate, fargs=(im, ax, Soln),
                         frames=np.arange(0, 20.e-6, 0.05e-6),
                         repeat=False)
    return anim


if __name__ == '__main__':
    animation = main()
    animation.save('Hexagon.mp4', fps=20)
    plt.close('all')
