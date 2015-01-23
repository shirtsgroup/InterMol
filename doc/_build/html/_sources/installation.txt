============
Installation
============

Install with pip
----------------

 $ pip install intermol

Dependencies
------------
To use InterMol, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Mac and CentOS machines. TravisCI is
        currently only set up to perform testing on Debian.

    `Python <http://python.org>`_ == 2.7
         Once our unit tests flesh out a bit more, we intend to add support
         for >=2.6.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        Numpy is the base package for numerical computing in python.

Testing your installation
-------------------------

InterMol uses `py.test` for unit testing. To run them simply type run the
following while in the base directory:

 $ py.test

We need a LOT more tests so any help here is especially welcome!
