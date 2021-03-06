# coding: utf-8

""" Utilities for Potential classes """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Third-party
import numpy as np

# Project
from .core import PotentialBase

__all__ = ['from_equation']

def _classnamify(s):
    s = [x.lower() for x in str(s).split()]
    words = []
    for word in s:
        words.append(word[0].upper() + word[1:])
    return "".join(words)

def from_equation(expr, vars, pars, name=None, hessian=False):
    r"""
    Create a potential class from an expression for the potential.

    .. note::

        This utility requires having `Sympy <http://www.sympy.org/>`_ installed.

    Parameters
    ----------
    expr : :class:`sympy.core.expr.Expr`, str
        Either a ``Sympy`` expression, or a string that can be converted to
        a ``Sympy`` expression.
    vars : iterable
        An iterable of variable names in the expression.
    pars : iterable
        An iterable of parameter names in the expression.
    name : str (optional)
        The name of the potential class returned.
    hessian : bool (optional)
        Generate a function to compute the Hessian.

    Returns
    -------
    CustomPotential : `~gary.potential.PotentialBase`
        A potential class that represents the input equation. To instantiate the
        potential, use just like a normal class with parameters.

    Examples
    --------
    Here we'll create a potential class for the harmonic oscillator
    potential, :math:`\Phi(x) = \frac{1}{2}\,k\,x^2`::

        >>> Potential = from_equation("1/2*k*x**2", vars="x", pars="k",
        ...                           name='HarmonicOscillator')
        >>> p1 = Potential(k=1.)
        >>> p1
        <HarmonicoscillatorPotential: k=1.00 (dimensionless)>

    The potential class (and object) is a fully-fledged subclass of
    `~gary.potential.PotentialBase` and therefore has many useful methods.
    For example, to integrate an orbit::

        >>> orbit = p1.integrate_orbit([1.,0], dt=0.01, nsteps=1000)

    """
    try:
        import sympy
        from sympy.utilities.lambdify import lambdify
    except ImportError:
        raise ImportError("sympy is required to use 'from_equation()' "
                          "potential class creation.")

    # convert all input to Sympy objects
    expr = sympy.sympify(expr)
    vars = [sympy.sympify(v) for v in vars]
    var_names = [v.name for v in vars]
    pars = [sympy.sympify(p) for p in pars]
    par_names = [p.name for p in pars]

    class MyPotential(PotentialBase):

        def __init__(self, units=None, **kwargs):
            self.parameters = kwargs
            for par in par_names:
                if par not in self.parameters:
                    raise ValueError("You must specify a value for "
                                     "parameter '{}'.".format(par))
            super(MyPotential,self).__init__(units)

    if name is not None:
        name = _classnamify(name)
        if "potential" not in name.lower():
            name = name + "Potential"
        MyPotential.__name__ = name

    # Energy / value
    valuefunc = lambdify(vars + pars, expr, dummify=False)
    def _value(self, w, t):
        kw = self.parameters.copy()
        for i,name in enumerate(var_names):
            kw[name] = w[i]
        return valuefunc(**kw)
    MyPotential._value = _value

    # Gradient
    gradfuncs = []
    for var in vars:
        gradfuncs.append(lambdify(vars + pars, sympy.diff(expr,var), dummify=False))

    def _gradient(self, w, t):
        kw = self.parameters.copy()
        for i,name in enumerate(var_names):
            kw[name] = w[i]
        return np.vstack([f(**kw) for f in gradfuncs])
    MyPotential._gradient = _gradient

    # Hessian
    if hessian:
        raise NotImplementedError("Hessian functionality doesn't exist yet...sorry.")
        hessfuncs = []

    return MyPotential
