# coding: utf-8

""" General dynamics utilities. """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
import astropy.coordinates as coord
from astropy import log as logger
import numpy as np

__all__ = ['Orbit']

class Orbit(object):
    """
    Represents the orbit of an object.

    Parameters
    ----------
    pos : array_like, :class:`~astropy.units.Quantity`
        Array of positions.
    vel : array_like, :class:`~astropy.units.Quantity`
        Array of velocities.
    t : array_like, :class:`~astropy.units.Quantity` (optional)
        Array of times.
    units : iterable (optional)
        Specify the unit system of the input orbit. This must be a tuple or iterable
        of :class:`~astropy.units.Unit` objects that define a unit system
    representation : :class:`~astropy.coordinates.BaseRepresentation` (optional)
        The representation of the input positions and velocities. By default,
        this is assumed to be :class:`~astropy.coordinates.CartesianRepresentation`
        so that the input coordinates are Cartesian, but this can be any of the
        valid Astropy representations.
    t_axis : int (optional)
        The index of the axis that represents the timesteps of the orbit.
    potential : :class:`~gary.potential.Potential` (optional)
        The potential that the orbit was integrated in.

    """

    def __init__(self, pos=None, vel=None, t=None, units=None,
                 representation=coord.CartesianRepresentation,
                 t_axis=0, potential=None):

        # TODO: figure out if tube or box (.orbit_type)
        pass

    def represent_as(self, representation):
        """
        Transform the representation or coordinate system of the orbit, for example,
        from Cartesian to Spherical.

        Parameters
        ----------
        representation : :class:`~astropy.coordinates.BaseRepresentation`
            The output representation.
        """
        pass


    def __getitem__(self, slyce):
        pass

    # TODO: for each representation shorthand, add an attribute to get out
    #       a new Orbit with that representation

    # TODO: should have attributes for whatever representation it's in,
    #       and also same with v in front for velocities

    def angular_momentum(self):
        pass

    def kinetic_energy(self):
        pass

    def potential_energy(self):
        pass

    def total_energy(self):
        pass

    def align_circulation_with_z(self):
        # TODO: returns copy
        pass

    def plot(self):
        pass
