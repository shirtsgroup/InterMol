#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import os
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

#####################################
VERSION = "0.2.0"
ISRELEASED = False
if ISRELEASED:
    __version__ = VERSION
else:
    __version__ = VERSION + '.dev0'
#####################################


with open('intermol/version.py', 'w') as version_file:
    version_file.write('version="{0}"\n'.format(__version__))

with open('__conda_version__.txt', 'w') as conda_version:
    conda_version.write(__version__)

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(['intermol'])
        sys.exit(errcode)

readme = open('README.md').read()

setup(
    name='intermol',
    version=__version__,
    description='InterMol is a conversion tool for molecular simulations.',
    long_description=readme,
    author='Christoph Klein, Christopher Lee, Ellen Zhong, and Michael Shirts',
    author_email='ctk3b@virginia.edu, ctl4f@virginia.edu, edz3fz@virginia.edu, '
                 'michael.shirts@virginia.edu',
    url='https://github.com/shirtsgroup/intermol',
    download_url='https://github.com/shirtsgroup/intermol/tarball/{}'.format(__version__),
    packages=find_packages(),
    package_dir={'intermol': 'intermol'},
    include_package_data=True,
    license="MIT",
    zip_safe=False,
    keywords='intermol',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
    test_suite='tests',
    cmdclass={'tests': PyTest,
    },
)
