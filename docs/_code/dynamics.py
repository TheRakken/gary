# coding: utf-8

""" ...explain... """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os, sys
import logging

# Third-party
import numpy as np
from astropy import log as logger
import astropy.units as u
import matplotlib.pyplot as plt

# Package
import streamteam.potential as sp
import streamteam.integrate as si
import streamteam.dynamics as sd

logger.setLevel(logging.DEBUG)

plot_path = os.path.abspath("../_static/dynamics")
if not os.path.exists(plot_path):
    os.mkdir(plot_path)

def main(overwrite=False):

    # only plot up to this index in the orbit plots
    nsteps = 500000
    plot_ix = nsteps//35

    # define an axisymmetric potential
    usys = (u.kpc, u.Msun, u.Myr)
    p = sp.LogarithmicPotential(v_c=0.15, r_h=0., phi=0.,
                                q1=1., q2=1., q3=0.85,  usys=usys)

    # initial conditions
    w0 = [8.,0.,0.,0.075,0.15,0.05]

    orbit_filename = os.path.join(plot_path, "orbits.npy")
    action_filename = os.path.join(plot_path, "actions.npy")

    if overwrite:
        if os.path.exists(orbit_filename):
            os.remove(orbit_filename)

        if os.path.exists(action_filename):
            os.remove(action_filename)

    if not os.path.exists(orbit_filename):
        # integrate an orbit in a axisymmetric potential
        acc = lambda t,x: p.acceleration(x)
        integrator = si.LeapfrogIntegrator(acc)
        t,w = integrator.run(w0, dt=1., nsteps=nsteps)

        # also integrate the orbit in the best-fitting isochrone potential
        m,b = sd.fit_isochrone(w, usys=usys)
        isochrone = sp.IsochronePotential(m=m, b=b, usys=usys)
        acc = lambda t,x: isochrone.acceleration(x)
        integrator = si.LeapfrogIntegrator(acc)
        iso_t,iso_w = integrator.run(w0, dt=1., nsteps=plot_ix)

        # cache the orbits
        np.save(orbit_filename, (t,w,iso_t,iso_w))
    else:
        t,w,iso_t,iso_w = np.load(orbit_filename)

    # plot a smaller section of the orbit in projections of XYZ
    fig = sd.plot_orbits(iso_w, marker=None, linestyle='-',
                         alpha=0.5, triangle=True, c='r')
    fig = sd.plot_orbits(w[:plot_ix], axes=fig.axes, marker=None, linestyle='-',
                         alpha=0.8, triangle=True, c='k')
    fig.savefig(os.path.join(plot_path, "orbit_xyz.png"))

    # plot a smaller section of the orbit in the meridional plane
    fig,ax = plt.subplots(1,1,figsize=(6,6))
    R = np.sqrt(w[:,0,0]**2 + w[:,0,1]**2)
    iso_R = np.sqrt(iso_w[:,0,0]**2 + iso_w[:,0,1]**2)
    ax.plot(iso_R, iso_w[:,0,2], marker=None, linestyle='-', alpha=0.5, c='r')
    ax.plot(R[:plot_ix], w[:plot_ix,0,2], marker=None, linestyle='-', alpha=0.8, c='k')
    ax.set_xlabel("$R$")
    ax.set_ylabel("$Z$", rotation='horizontal')
    fig.savefig(os.path.join(plot_path, "orbit_Rz.png"))

    if not os.path.exists(action_filename):
        # compute the actions and angles for the orbit
        actions,angles,freqs = sd.cross_validate_actions(t, w[:,0], N_max=6, nbins=100, usys=usys)

        # now compute for the full time series
        r = sd.find_actions(t, w[:,0], N_max=6, usys=usys, return_Sn=True)
        full_actions,full_angles,full_freqs = r[:3]
        Sn,dSn_dJ,nvecs = r[3:]

        np.save(action_filename, (actions,angles,freqs) + r)
    else:
        r = np.load(action_filename)
        actions,angles,freqs = r[:3]
        full_actions,full_angles,full_freqs = r[3:6]
        Sn,dSn_dJ,nvecs = r[6:]

    # deviation of actions computed from subsample to median value
    fig,axes = plt.subplots(1,3,figsize=(12,5),sharey=True,sharex=True)
    bins = np.linspace(-0.1,0.1,20)
    for i in range(3):
        axes[i].set_title("$J_{}$".format(i+1), y=1.02)
        axes[i].plot((actions[:,i] - np.median(actions[:,i])) / np.median(actions[:,i])*100.,
                     marker='.', linestyle='none')
        axes[i].set_ylim(-.11,.11)

    axes[0].set_yticks((-0.1,-0.05,0.,0.05,0.1))
    axes[0].set_yticklabels(["{}%".format(tck) for tck in axes[0].get_yticks()])

    axes[1].set_xlabel("subsample index")

    fig.suptitle("Deviation from median action value", fontsize=20)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_path, "action_hist.png"))

    # deviation of frequencies computed from subsample to median value
    fig,axes = plt.subplots(1,3,figsize=(12,5),sharey=True,sharex=True)
    bins = np.linspace(-0.1,0.1,20)
    for i in range(3):
        axes[i].set_title(r"$\Omega_{}$".format(i+1), y=1.02)
        axes[i].plot((freqs[:,i] - np.median(freqs[:,i])) / np.median(freqs[:,i])*100.,
                     marker='.', linestyle='none')
        axes[i].set_ylim(-.11,.11)

    axes[0].set_yticks((-0.1,-0.05,0.,0.05,0.1))
    axes[0].set_yticklabels(["{}%".format(tck) for tck in axes[0].get_yticks()])

    axes[1].set_xlabel("subsample index")

    fig.suptitle("Deviation from median frequency value", fontsize=20)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_path, "freq_hist.png"))

    # --------------------------------------------------------
    # now going to plot toy actions and solved actions

    # fit isochrone potential
    m,b = sd.fit_isochrone(w, usys=usys)
    isochrone = sp.IsochronePotential(m=m, b=b, usys=usys)
    actions,angles = isochrone.action_angle(w[:,0,:3],w[:,0,3:])

    fig,axes = plt.subplots(1,3,figsize=(12,5),sharey=True,sharex=True)
    for i in range(3):
        computed_action = full_actions[i]
        axes[i].plot(t/1000., (actions[:,i]-computed_action)/computed_action*100,
                     marker=None, alpha=0.5, label='toy action', lw=1.5)
        axes[i].axhline(0., lw=1., zorder=-1, c='#31a354')
        axes[i].set_title("$J_{}$".format(i+1), y=1.02)

    axes[1].set_xlabel("time [Gyr]")

    fig.suptitle("Percent deviation from estimated action", fontsize=20)
    axes[0].legend(fontsize=16)
    axes[0].set_yticks((-50,-25,0,25,50))
    axes[0].set_yticklabels(["{}%".format(tck) for tck in axes[0].get_yticks()])

    dt = t[1]-t[0]
    axes[0].set_xlim(0.,3.)
    axes[0].set_ylim(-52,52)
    fig.tight_layout()
    fig.savefig(os.path.join(plot_path,"toy_computed_actions.png"))

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Create logger
    logger = logging.getLogger(__name__)
    ch = logging.StreamHandler()
    formatter = logging.Formatter("%(name)s / %(levelname)s / %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Define parser object
    parser = ArgumentParser(description="")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help="Be chatty! (default = False)")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
                        default=False, help="Be quiet! (default = False)")

    parser.add_argument("-o", dest="overwrite", action="store_true", default=False,
                        help="Overwrite generated files.")

    args = parser.parse_args()

    # Set logger level based on verbose flags
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    elif args.quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    main(overwrite=args.overwrite)