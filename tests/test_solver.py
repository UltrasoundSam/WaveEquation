#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import numpy as np
import pytest

from waveequation import WaveSolver


def test_velocity_size():
    # Set up domain
    x = np.arange(100)
    y = np.arange(150)
    correct_size = (len(x), len(y))
    A = np.ones(correct_size)
    t0 = np.zeros(correct_size)

    # Make velocity array incorrect size
    c = np.ones([100, 130])*3230

    # Check error is caught
    with pytest.raises(ValueError):
        WaveSolver(x, y, c, dt=2.5e-9, params=(A, 2e6, 4e-7, t0))
