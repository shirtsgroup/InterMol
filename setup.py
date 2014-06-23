#! /usr/bin/python

from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension, find_packages

if __name__ == '__main__':
    RELEASE = "0.1-dev" 
    LONG_DESCRIPTION = ""
    CLASSIFIERS = ['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ]

    extensions = []

    setup(name              = 'intermol',
          version           = RELEASE,
          description       = 'Conversion tool for molecular simulations',
          author            = 'Christoph Klein, Christopher Lee, Ellen Zhong, and Michael Shirts',
          author_email      = '<ctk3b@virginia.edu>, <ctl4f@virginia.edu>, <edz3fz@virginia.edu>, <michael.shirts@virginia.edu>',
          url               = 'https://github.com/shirtsgroup/InterMol',
          requires          = [],
          provides          = ['intermol'],
          packages          = [ 'intermol',
                                'intermol.desmond_extension',
                                'intermol.gromacs_extension',
                                'intermol.forces',
                                'intermol.lammps_extension',
                                'intermol.types',
                                'intermol.unit',],
          package_dir       = {'intermol': 'intermol'},
          ext_package       = 'intermol',
          license           = 'xxx',
          ext_modules       = extensions,
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          #cmdclass          = cmdclass,
          # all standard requirements are available through PyPi and
          # typically can be nstalled without difficulties through setuptools
          install_requires = [],
          # extras can be difficult to install through setuptools and/or
          # you might prefer to use the version available through your
          # packaging system
          extras_require = {                },
          test_suite = "InterMolTests",
          tests_require = ['nose>=0.10',
                           ],
          zip_safe = False,     # as a zipped egg the *.so files are not found (at least in Ubuntu/Linux)
          )
