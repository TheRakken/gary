# coding: utf-8

""" Test core dynamics.  """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os
import logging

# Third-party
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import log as logger

# Project
from ..core import *
from ..plot import plot_orbits
from ..orbit import Orbit
from ...potential import LogarithmicPotential
from ...units import galactic

logger.setLevel(logging.DEBUG)

plot_path = "plots/tests/dynamics"
if not os.path.exists(plot_path):
    os.makedirs(plot_path)

def test_angular_momentum():
    # TODO: need to try different representations...

    qs = [[1.,0.,0.],[1.,0.,0.],[0.,1.,0.]]
    ps = [[0.,0.,1.],[0.,1.,0.],[0.,0.,1.]]
    Ls = [[0., -1, 0],[0.,0.,1.],[1.,0.,0.]]

    for q,p,L in zip(qs,ps,Ls):
        L1 = angular_momentum(pos=q*u.dimensionless_unscaled,
                              vel=p*u.dimensionless_unscaled)
        np.testing.assert_allclose(L1, L)

        L2 = angular_momentum(pos=q*u.kpc,
                              vel=p*u.km/u.s)
        np.testing.assert_allclose(L2, L)

        orbit = Orbit(pos=q*u.kpc, vel=p*u.km/u.s)
        L3 = angular_momentum(orbit)
        np.testing.assert_allclose(L3, L)

    # now try with arrays
    q,p,L = map(np.transpose, map(np.array, (qs,ps,Ls)))
    L1 = angular_momentum(pos=q*u.dimensionless_unscaled,
                          vel=p*u.dimensionless_unscaled)

    np.testing.assert_allclose(L1, L)

    L2 = angular_momentum(pos=q*u.kpc,
                          vel=p*u.km/u.s)
    np.testing.assert_allclose(L2, L)

    orbit = Orbit(pos=q*u.kpc, vel=p*u.km/u.s)
    L3 = angular_momentum(orbit)
    np.testing.assert_allclose(L3, L)

    # try with multiple orbits
    orbit = Orbit(pos=np.random.random(size=(3,1000,5))*u.kpc,
                  vel=np.random.random(size=(3,1000,5))*u.km/u.s)
    L1 = angular_momentum(orbit)
    assert L1.shape == (3,1000,5)

# ----------------------------------------------------------------------------

def test_classify_orbit():
    # TODO: need to try different representations...
    t = np.linspace(0, 10., 100)

    # loopy orbit
    xyz = np.vstack((np.cos(t), np.sin(t), t*0.)) * u.au
    vxyz = np.vstack((-np.sin(t), np.cos(t), t*0.)) * u.au/u.yr
    loop = classify_orbit(pos=xyz, vel=vxyz)
    assert loop.sum() == 1
    assert loop.shape == (3,)

    # two loopy orbits
    xyz2 = np.vstack((np.cos(2*t), t*0., np.sin(2*t))) * u.au
    vxyz2 = np.vstack((-2*np.sin(2*t), t*0., 2*np.cos(2*t))) * u.au/u.yr

    pos = np.zeros((3,100,2))
    pos[...,0] = xyz.value
    pos[...,1] = xyz2.value
    pos = pos*xyz.unit

    vel = np.zeros((3,100,2))
    vel[...,0] = vxyz.value
    vel[...,1] = vxyz2.value
    vel = vel*vxyz.unit
    loop = classify_orbit(pos=pos, vel=vel)
    np.testing.assert_allclose(loop, np.array([[0,0],[0,1],[1,0]]))

    # box-like orbit
    xyz = np.vstack((np.cos(t), np.cos(2*t + 0.5), np.cos(4*t + 1.5))) * u.au
    vxyz = -np.vstack((np.sin(t), 2*np.sin(2*t + 0.5), 4*np.sin(4*t + 1.5))) * u.au/u.yr
    loop = classify_orbit(pos=xyz, vel=vxyz)

    assert loop.sum() == 0
    assert loop.shape == (3,)

# ----------------------------------------------------------------------------

def test_align_circulation_single():

    potential = LogarithmicPotential(v_c=1., r_h=0.14, q1=1., q2=0.9, q3=1.,
                                     units=galactic)
    w0 = np.array([[0.,1.,0.,0.,0.,0.5],  # loop around x axis
                   [1.,0.,0.,0.,0.,0.5],  # loop around y axis
                   [1.,0.,0.,0.,0.5,0.],  # loop around z axis
                   [0.8,0.4,0.,0.,0.1,0.]])  # box

    t,w = potential.integrate_orbit(w0, dt=0.05, nsteps=10000)

    for i in range(w.shape[1]):
        circ = classify_orbit(w[:,i])
        new_w = align_circulation_with_z(w[:,i], circ)
        new_circ = classify_orbit(new_w)

        if i == 3:
            assert np.sum(new_circ) == 0
        else:
            assert new_circ[2] == 1.

def test_align_circulation_many():

    potential = LogarithmicPotential(v_c=1., r_h=0.14, q1=1., q2=0.9, q3=1.,
                                     units=galactic)
    w0 = np.array([[0.,1.,0.,0.,0.,0.5],  # loop around x axis
                   [1.,0.,0.,0.,0.,0.5],  # loop around y axis
                   [1.,0.,0.,0.,0.5,0.],  # loop around z axis
                   [0.8,0.4,0.,0.,0.1,0.]])  # box
    names = ['xloop', 'yloop', 'zloop', 'box']

    t,w = potential.integrate_orbit(w0, dt=0.05, nsteps=10000)
    fig = plot_orbits(w, linestyle='none', alpha=0.1)
    fig.savefig(os.path.join(plot_path, "align_circulation_orbits_init.png"))

    circ = classify_orbit(w)
    assert circ.shape == (4,3)

    new_w = align_circulation_with_z(w, circ)
    fig = plot_orbits(new_w, linestyle='none', alpha=0.1)
    fig.savefig(os.path.join(plot_path, "align_circulation_orbits_post.png"))

    new_circ = classify_orbit(new_w)
    assert np.all(new_circ[:3,2] == 1.)
