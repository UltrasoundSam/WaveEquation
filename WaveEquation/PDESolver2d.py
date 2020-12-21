#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:07:45 2020

@author: sam
"""
import numpy as np

class WaveSolver:
    '''
    Class for numerically solving the wave equation
    '''
    def __init__(self, x, y, c, u_init=None, dt=None, params=None):
        '''
        Initialise the problem by creating grid, etc.
        Inputs:
            x       -   Array describing x-coordinates of grid
            y       -   Array describing y-coordinates of grid
            c       -   2D array describing wave velocity at each x,y location
            u_init  -   2D grid describing the initial condition
            dt      -   User defined timestep, must be less than the minimum
                        time it takes for the wave to cross a meshpoint.
            params  -   Tuple containing parameters used to control excitation
                        function; specfically (A, freq, width, timedelay):
                            A       - Amplitude [array same size as c]
                            freq    - Frequency of wave
                            width   - Width of Gaussian window
                            t0      - Time delay on function [int or array]
        '''
        # Read in information about domain, and create grid
        self.x = x
        self.dx = x[2] - x[1]
        self.y = y
        self.dy = y[2] - y[1]

        # Create spatial grid - the current, and two previous time steps
        self.u = np.zeros([len(x), len(y)])
        if not u_init:
            self.u_1 = np.zeros(self.u.shape)
        else:
            self.u_1 = u_init
        self.u_2 = np.zeros(self.u.shape)

        if not params:
            self.set_functionvalues(1, 2e6, 4e-7, 2e-6)
        else:
            self.set_functionvalues(*params)

        # Any active driving function
        self.f = np.zeros(self.u.shape)
        self.f[self.mask] = np.real(self.gaussian(time=0,
                                                  A=self.amp[self.mask],
                                                  freq=self.frequency,
                                                  sigma=self.sigma,
                                                  t0=self.t0[self.mask]))

        # Calculate the timestep
        self.cmax = c[c != 0].max()
        if not dt:
            self.dt = 0.8 * min(self.dx, self.dy)/self.cmax
        else:
            assert dt < (min(self.dx, self.dy)/self.cmax), "Timestep too large."
            self.dt = dt

        # Calculate the Courants
        self.courantx = (c*self.dt/self.dx)**2
        self.couranty = (c*self.dt/self.dy)**2

        # Calculate the first timestep
        xgradient = self.u_1[:-2, 1:-1] - 2*self.u_1[1:-1, 1:-1] + self.u_1[2:, 1:-1]
        ygradient = self.u_1[1:-1, :-2] - 2*self.u_1[1:-1, 1:-1] + self.u_1[1:-1, 2:]
        self.u[1:-1, 1:-1] = self.u_1[1:-1, 1:-1] \
                            + 0.5*(self.courantx[1:-1, 1:-1]*xgradient + self.couranty[1:-1, 1:-1]*ygradient) \
                            + self.f[1:-1, 1:-1]

        # Copy data over for next timestep calculation
        self.u_2 = self.u_1.copy()
        self.u_1 = self.u.copy()
        self.t = self.dt

    def set_timestep(self, value):
        '''
        Sets timestep to given value, assuming it meets stability
        criteria (that is, it is less than the time it takes a wave to travel
        from one element to another)
        '''
        assert value < (min(self.dx, self.dy)/self.cmax), "Timestep too large."
        self.dt = value

    def set_functionvalues(self, amplitude, frequency, width, timedelay):
        '''
        Allows the user to set certain parameters that can be used to control
        the excitation function. It also creates a mask that determines the
        active parts of the mesh to reduce calculations.
        Inputs:
            amplitude   -   Amplitude of wave
            frequency   -   Frequency of wave
            width       -   Width of wavepacket
            timedelay   -   Time delay of pulse
        '''
        # Determining the active points in the array
        self.mask = (amplitude != 0)

        # Assigning parameters
        self.amp = amplitude
        self.frequency = frequency
        self.sigma = width
        # Making sure that timedelay has the same shape as amplitude
        try:
            timedelay[self.mask]
        except TypeError:
            # Must be just an integer, so make it the right shape
            timebuff = timedelay
            timedelay = np.zeros(self.mask.shape)
            timedelay[self.mask] = timebuff
        self.t0 = timedelay

    def gaussian(self, time, A, freq, sigma, t0):
        '''
        Function that creates gaussian windowed sine function, with the
        following parameters:
            time    -   Time at which function is evaluated
            A       -   Maximum amplitude of function
            freq    -   Frequency
            Sigma   -   Width of wavepacket
            t0      -   Time delay
        '''
        return A*np.exp(-(time-t0)**2/sigma**2)*np.exp(2j*np.pi*freq*(time-t0))

    def solve_step(self):
        '''
        Iterates forward by one step
        '''
        self.f[self.mask] = np.real(self.gaussian(time=self.t,
                                                  A=self.amp[self.mask],
                                                  freq=self.frequency,
                                                  sigma=self.sigma,
                                                  t0=self.t0[self.mask]))

        # Calculate gradients
        xgradient = self.u_1[:-2, 1:-1] - 2*self.u_1[1:-1, 1:-1] + self.u_1[2:, 1:-1]
        ygradient = self.u_1[1:-1, :-2] - 2*self.u_1[1:-1, 1:-1] + self.u_1[1:-1, 2:]

        # Calculate next step
        self.u[1:-1, 1:-1] = 2*self.u_1[1:-1, 1:-1] - self.u_2[1:-1, 1:-1] \
                            + self.courantx[1:-1, 1:-1]*xgradient \
                            + self.couranty[1:-1, 1:-1]*ygradient \
                            + self.f[1:-1, 1:-1]
        # Update previous grids
        self.u_2 = self.u_1.copy()
        self.u_1 = self.u.copy()
        self.t += self.dt

    def get_snapshot(self, decimate=1):
        '''
        Extract current wavefield. Can also extract every nth value in x and
        y by setting the decimate parameter to n.
        '''
        return self.u[::decimate, ::decimate]
