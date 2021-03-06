.. include:: references.txt

.. _gary-install:

************
Installation
************

Dependencies
============

This packages has the following dependencies:

- `Python`_ >= 2.7

- `Numpy`_ >= 1.8

- `Cython <http://www.cython.org/>`_: >= 0.21

- `Astropy`_ >= 1.0

- `PyYAML`_ >= 3.10

You can use ``pip`` or ``conda`` to install these automatically.

Optional Dependencies
---------------------

- PyGaia (`pip install pygaia`)

- `Sympy`_

- `emcee`_ (only for running some tests)

Installing
==========

Development repository
----------------------

The latest development version of gary can be cloned from
`GitHub <https://github.com/>`_ using ``git``::

   git clone git://github.com/adrn/gary.git

Building and Installing
-----------------------

To build the project (from the root of the source tree, e.g., inside
the cloned ``gary`` directory)::

    python setup.py build

To install the project::

    python setup.py install
