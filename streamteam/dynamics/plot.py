# coding: utf-8

""" ...explain... """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os, sys

# Third-party
import matplotlib.pyplot as plt
import numpy as np

__all__ = ['plot_orbits']

def plot_orbits(w, ix=None, axes=None, triangle=False, **kwargs):
    """
    TODO:
    """

    if triangle and axes is None:
        fig,axes = plt.subplots(2,2,figsize=(12,12),sharex='col',sharey='row')
        axes[0,1].set_visible(False)
        axes = axes.flat
        axes = [axes[0],axes[2],axes[3]]

    elif triangle and axes is not None:
        try:
            axes = axes.flat
        except:
            pass

        if len(axes) == 4:
            axes = [axes[0],axes[2],axes[3]]

    elif not triangle and axes is None:
        fig,axes = plt.subplots(1,3,figsize=(12,5),sharex=True,sharey=True)

    if ix is not None:
        ixs = [ix]
    else:
        ixs = range(w.shape[1])

    for ii in ixs:
        axes[0].plot(w[:,ii,0], w[:,ii,1], **kwargs)
        axes[1].plot(w[:,ii,0], w[:,ii,2], **kwargs)
        axes[2].plot(w[:,ii,1], w[:,ii,2], **kwargs)

    if triangle:
        # HACK: until matplotlib 1.4 comes out, need this
        axes[0].set_ylim(axes[0].get_xlim())
        axes[2].set_xlim(axes[0].get_ylim())

        axes[0].set_ylabel("Y")
        axes[1].set_xlabel("X")
        axes[1].set_ylabel("Z")
        axes[2].set_xlabel("Y")

    else:
        axes[0].set_xlabel("X")
        axes[0].set_ylabel("Y")

        axes[1].set_xlabel("X")
        axes[1].set_ylabel("Z")

        axes[2].set_xlabel("Y")
        axes[2].set_ylabel("Z")

    if not triangle:
        axes[0].figure.tight_layout()

    return axes[0].figure