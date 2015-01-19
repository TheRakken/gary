# coding: utf-8

""" General dynamics utilities. """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import astropy.units as u
import astropy.coordinates as coord
from astropy import log as logger
import numpy as np

# Project
from ..coordinates import velocity_transforms as vtrans

__all__ = ['Orbit']

def readonly_property(obj, name, value):
    setattr(obj.__class__, name, property(lambda obj: obj.__dict__["_" + name]))
    setattr(obj, "_" + name, value)

class Orbit(object):
    """
    Represents the orbit of an object.

    Parameters
    ----------
    pos : array_like, :class:`~astropy.units.Quantity`
        Array of positions. axis=0 is assumed to be the coordinate axis and should have
        length=3. axis=1 is assumed to be the time axis, but is optional (for initial
        conditions). Any further axes should be handled properly.
    vel : array_like, :class:`~astropy.units.Quantity`
        Array of velocities. axis=0 is assumed to be the coordinate axis and should have
        length=3. axis=1 is assumed to be the time axis, but is optional (for initial
        conditions). Any further axes should be handled properly.
    t : array_like, :class:`~astropy.units.Quantity` (optional)
        Array of times.
    unitsys : iterable (optional)
        Specify the unit system of the input orbit. This must be a tuple or iterable
        of :class:`~astropy.units.Unit` objects that define a unit system
    representation : :class:`~astropy.coordinates.BaseRepresentation` (optional)
        The representation of the input positions and velocities. By default,
        this is assumed to be :class:`~astropy.coordinates.CartesianRepresentation`
        so that the input coordinates are Cartesian, but this can be any of the
        valid Astropy representations.
    potential : :class:`~gary.potential.Potential` (optional)
        The potential that the orbit was integrated in.

    """

    def __init__(self, pos, vel, t=None, unitsys=None,
                 Representation=coord.CartesianRepresentation,
                 t_axis=1, potential=None):

        if not hasattr(pos, "unit") or pos.shape[0] != 3:
            raise TypeError("Input position must be an Astropy Quantity object "
                            "with an axis=0 of length 3 (e.g., x,y,z or rho,phi,z).")

        if not hasattr(vel, "unit") or vel.shape[0] != 3:
            raise TypeError("Input velocity must be an Astropy Quantity object "
                            "with an axis=0 of length 3 (e.g., vx,vy,vz).")

        if pos.shape != vel.shape:
            raise ValueError("Position and velocity must have the same shape "
                             " ({} vs. {})".format(pos.shape, vel.shape))

        if t is not None and not hasattr(t, "unit") or t.size != pos.shape[1]:
            raise TypeError("Input time must be an Astropy Quantity object and match "
                            "the size of the time axis in position and velocity.")

        self.pos = pos
        self.vel = vel
        self.t = t
        self.unitsys = unitsys
        self.Representation = Representation
        self.potential = potential

        for i,name in enumerate(Representation.attr_classes.keys()):
            readonly_property(self, name, self.pos[i])
            readonly_property(self, "v"+name, self.vel[i])

    def represent_as(self, Representation):
        """
        Transform the representation or coordinate system of the orbit, for example,
        from Cartesian to Spherical.

        Parameters
        ----------
        Representation : :class:`~astropy.coordinates.BaseRepresentation`
            The output representation class.
        """

        # first transform the position
        new_pos = self.Representation(self.pos).represent_as(Representation)

        # now find the function to transform the velocity
        _func_name = "{}_to_{}".format(self.Representation.get_name(),
                                       Representation.get_name())
        v_func = getattr(vtrans, _func_name)
        new_vel = v_func(self.pos, self.vel)

        return Orbit(pos=new_pos, vel=new_vel, t=self.t, unitsys=self.unitsys,
                     Representation=Representation, potential=self.potential)

    def __repr__(self):
        return "<Orbit"

    def __str__(self):
        pass

    def __getitem__(self, slyce):
        pass

    # TODO: for each representation shorthand, add an attribute to get out
    #       a new Orbit with that representation

    # TODO: should have attributes for whatever representation it's in,
    #       and also same with v in front for velocities

    @property
    def orbit_type(self):
        # TODO: figure out if tube or box (.orbit_type)
        pass

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
