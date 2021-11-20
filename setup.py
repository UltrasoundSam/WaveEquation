#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup
import WaveEquation
import os

"""Utility function to read the README file.
Used for the long_description.  It's nice, because now 1) we have a top level
README file and 2) it's easier to type in the README file than to put a raw
string in below ..."""

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='WaveEquation',
    version='0.1.0',
    description='A simple finite difference solver for the 2D wave-equation',
    author='Sam Hill',
    author_email='sam.j.hill@protonmail.com',
    packages=['WaveEquation'],
    install_requires = ['numpy', 'matplotlib']
     )
