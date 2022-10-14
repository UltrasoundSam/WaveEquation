# WaveEquation
Simple script that numerically solves the wave equation using finite difference methods. Nothing extravagent, just something I knocked up quickly over the Christmas break as I was a little bored. 

The solver takes any 2D grid of wave speed values, and solve the 2D linear partial differential equation in order to describe how a wave travels through the domain given:

- Initial condition of the wavefield
- Any driving forces that may exist anywhere in the domain
- Fixed (Dirichlet) boundary conditions

# Example
As an example, see the output from a simulation of an angled beam (25 degrees) travelling through a medium with a speed of 3000 m/s. The top-right quadrant comprises of a material with a wavespeed of 1500 m/s. Diffraction, refraction and reflection can all be seen. 

![Wave Example](img/demo.gif)

![Tests](https://github.com/UltrasoundSam/WaveEquation/actions/workflows/tests.yml/badge.svg)


