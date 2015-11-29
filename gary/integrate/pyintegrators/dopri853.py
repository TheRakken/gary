# coding: utf-8

""" Wrapper around SciPy DOPRI853 integrator. """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
from scipy.integrate import ode

# Project
from ..core import Integrator
from ..timespec import _parse_time_specification

__all__ = ["DOPRI853Integrator"]

class DOPRI853Integrator(Integrator):
    r"""
    This provides a wrapper around `Scipy`'s implementation of the
    Dormand-Prince 85(3) integration scheme.

    .. seealso::

        - Numerical recipes (Dopr853)
        - http://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method

    Parameters
    ----------
    func : func
        A callable object that computes the phase-space coordinate
        derivatives with respect to the independent variable at a point
        in phase space.
    func_args : tuple (optional)
        Any extra arguments for the function.

    """

    def __init__(self, func, func_args=(), **kwargs):
        super(DOPRI853Integrator, self).__init__(func, func_args)
        self._ode_kwargs = kwargs

    def run(self, w0, mmap=None, **time_spec):
        """
        Run the integrator starting at the given coordinates and momenta
        (or velocities) and a time specification. The initial conditions
        `w0` should have shape `(nparticles, ndim)` or `(ndim,)` for a
        single orbit.

        There are a few combinations of keyword arguments accepted for
        specifying the timestepping. For example, you can specify a fixed
        timestep (`dt`) and a number of steps (`nsteps`), or an array of
        times. See **Other Parameters** below for more information.

        Parameters
        ==========
        w0 : array_like
            Initial conditions.

        Other Parameters
        ================
        dt, nsteps[, t1] : (numeric, int[, numeric])
            A fixed timestep dt and a number of steps to run for.
        dt, t1, t2 : (numeric, numeric, numeric)
            A fixed timestep dt, an initial time, and a final time.
        t : array_like
            An array of times to solve on.

        Returns
        =======
        times : array_like
            An array of times.
        w : array_like
            The array of positions and momenta (velocities) at each time in
            the time array. This array has shape `(Ntimes,Norbits,Ndim)`.

        """

        # generate the array of times
        times = _parse_time_specification(**time_spec)
        nsteps = len(times)-1

        w0, ws = self._prepare_ws(w0, mmap, nsteps)
        nparticles, ndim = w0.shape

        # need this to do resizing, and to handle func_args because there is some
        #   issue with the args stuff in scipy...
        def func_wrapper(t,x):
            _x = x.reshape((nparticles,ndim))
            return self.F(t, _x, *self._func_args).reshape((nparticles*ndim,))

        self._ode = ode(func_wrapper, jac=None)
        self._ode = self._ode.set_integrator('dop853', **self._ode_kwargs)

        # create the return arrays
        ws[:,0] = w0

        # make 1D
        w0 = w0.reshape((nparticles*ndim,))

        # set the initial conditions
        self._ode.set_initial_value(w0, times[0])

        # Integrate the ODE(s) across each delta_t timestep
        k = 1
        while self._ode.successful() and k < (nsteps+1):
            self._ode.integrate(times[k])
            outy = self._ode.y
            ws[:,k] = outy.reshape(nparticles,ndim)
            k += 1

        if not self._ode.successful():
            raise RuntimeError("ODE integration failed!")

        if ws.shape[-1] == 1:
            ws = ws[...,0]
        return times, ws
