# coding: utf-8

""" General dynamics utilities. """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import numpy as np
from astropy import log as logger
import astropy.coordinates as coord

# Project
from .orbit import Orbit

__all__ = ['angular_momentum', 'classify_orbit', 'align_circulation_with_z']

def angular_momentum(orbit=None, **kwargs):
    r"""
    Compute the angular momentum vector(s) of a single or array of positions
    and conjugate momenta, ``q`` and ``p``.

    .. math::

        \boldsymbol{L} = \boldsymbol{q} \times \boldsymbol{p}

    Parameters
    ----------
    orbit : :class:`~gary.dynamics.Orbit`
        Orbit instance.
    **keyword_arguments
        Instead of passing in an `~gary.dynamics.Orbit`, it is also possible to pass in
        a set of keyword arguments that fully specifies an `~gary.dynamics.Orbit` instance.
        For example, you may pass in ``pos=``, ``vel=``, and (optionally) a representation
        class. Any of the arguments accepted by the `~gary.dynamics.Orbit` initializer are
        also accepted here.

    Returns
    -------
    L : :class:`~astropy.units.Quantity`
        Array of angular momentum vectors.

    Examples
    --------

        >>> import astropy.units as u
        >>> import gary.dynamics as gdyn
        >>> orb = gdyn.Orbit(pos=[1., 0, 0]*u.au, vel=[0, 2*np.pi, 0]*u.au/u.yr)
        >>> gdyn.angular_momentum(orb)
        <Quantity [ 0.        , 0.        , 6.28318531] AU2 / yr>

    """

    if orbit is None:
        orbit = Orbit(**kwargs)

    cart_orbit = orbit.represent_as(coord.CartesianRepresentation)
    L = np.cross(cart_orbit.pos, cart_orbit.vel,
                 axisa=0, axisb=0, axisc=0) * cart_orbit.pos.unit * cart_orbit.vel.unit

    return L

def classify_orbit(orbit=None, **kwargs):
    """
    Determine whether an orbit or series of orbits is a Box or Tube orbit by
    figuring out whether there is a change of sign of the angular momentum
    about an axis. Returns a 2D array with 3 integers per orbit point such that:

    - Box and boxlet = [0,0,0]
    - z-axis (short-axis) tube = [0,0,1]
    - x-axis (long-axis) tube = [1,0,0]

    Parameters
    ----------
    orbit : :class:`~gary.dynamics.Orbit`
        Orbit instance.
    **keyword_arguments
        Instead of passing in an `~gary.dynamics.Orbit`, it is also possible to pass in
        a set of keyword arguments that fully specifies an `~gary.dynamics.Orbit` instance.
        For example, you may pass in ``pos=``, ``vel=``, and (optionally) a representation
        class. Any of the arguments accepted by the `~gary.dynamics.Orbit` initializer are
        also accepted here.

    Returns
    -------
    circulation : :class:`numpy.ndarray`
        An array that specifies whether there is circulation about any of
        the axes of the input orbit. For a single orbit, will return a
        1D array, but for multiple orbits, the shape will be (len(w), 3).

    Examples
    --------

        >>> import astropy.units as u
        >>> import gary.dynamics as gdyn
        >>> t = np.linspace(0., 10., 100)
        >>> xyz = np.vstack((np.cos(t), np.sin(t), t*0.)) * u.au
        >>> vxyz = np.vstack((-np.sin(t), np.cos(t), t*0.)) * u.au/u.yr
        >>> orb = gdyn.Orbit(pos=xyz, vel=vxyz)
        >>> gdyn.classify_orbit(orb)
        np.array([0,0,1])


    """

    if orbit is None:
        orbit = Orbit(**kwargs)

    # get angular momenta
    Ls = angular_momentum(orbit)

    # initial angular momentum
    L0 = Ls[:,0]

    # see if at any timestep the sign has changed
    if len(orbit.shape) > 1:
        loop = np.ones((3,) + orbit.shape[1:])
        single_orbit = False
    else:
        loop = np.ones((3,1))
        single_orbit = True

    for ii in range(3):
        cnd = (np.sign(L0[ii]) != np.sign(Ls[ii,1:])) | \
              (np.abs(Ls[ii,1:].value) < 1E-13)

        ix = np.atleast_1d(np.any(cnd, axis=0))
        loop[ii][ix] = 0

    loop = loop.astype(int)
    if single_orbit:
        return loop.reshape((3,))
    else:
        return loop

def align_circulation_with_z(orbit=None, circulation=None, **kwargs):
    """
    If the input orbit is a tube orbit, this function aligns the circulation
    axis with the z axis.

    Parameters
    ----------
    orbit : :class:`~gary.dynamics.Orbit`
        Orbit instance.
    circulation : array_like (optional)
        Array of bits that specify the axis about which the orbit circulates.
        See the documentation for `~gary.dynamics.classify_orbit`. If not
        provided, will run this function on the input orbit(s).
    **keyword_arguments
        Instead of passing in an `~gary.dynamics.Orbit`, it is also possible to pass in
        a set of keyword arguments that fully specifies an `~gary.dynamics.Orbit` instance.
        For example, you may pass in ``pos=``, ``vel=``, and (optionally) a representation
        class. Any of the arguments accepted by the `~gary.dynamics.Orbit` initializer are
        also accepted here.

    Returns
    -------
    new_orbit : :class:`~gary.dynamics.Orbit`
        A copy of the input `~gary.dynamics.Orbit` with the circulation
        aligned with the z axis. The representation of the output orbit
        will always be Cartesian.
    """

    if orbit is None:
        orbit = Orbit(**kwargs)

    if circulation is None:
        circulation = classify_orbit(orbit)

    orbit = orbit.represent_as(coord.CartesianRepresentation)

    if orbit.norbits == 1:
        issingle = True
        orbit = orbit[...,np.newaxis]
    else:
        issingle = False

    new_pos = orbit.pos.copy()
    new_vel = orbit.vel.copy()
    for ix in range(orbit.norbits):
        if circulation[2,ix] == 1 or np.all(circulation[:,ix] == 0):
            # already circulating about z or box orbit
            continue

        if sum(circulation[:,ix]) > 1:
            raise ValueError("Circulation about x and y axes - are you sure the "
                             "orbit has been integrated for long enough?")

        if circulation[0,ix] == 1:
            circ = 0
        elif circulation[1,ix] == 1:
            circ = 1
        else:
            raise RuntimeError("Should never get here...")

        new_pos[circ,:,ix] = orbit.pos[2,:,ix]
        new_pos[2,:,ix] = orbit.pos[circ,:,ix]
        new_vel[circ,:,ix] = orbit.vel[2,:,ix]
        new_vel[2,:,ix] = orbit.vel[circ,:,ix]

    new_orbit = orbit.copy()
    new_orbit.pos = new_pos
    new_orbit.vel = new_vel

    if issingle:
        new_orbit = new_orbit[...,0]

    return new_orbit
