#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import intermol.version


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(['intermol', '-s'])
        sys.exit(errcode)

with open('requirements.txt') as reqs:
    requirements_lines = [line.strip() for line in reqs]
reqs = list(filter(None, requirements_lines))

readme = open('README.md').read()

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

setup(
    name='intermol',
    version=intermol.version.short_version,
    description='InterMol is a conversion tool for molecular simulations.',
    long_description=readme,
    author='Christoph Klein, Christopher Lee, Ellen Zhong, and Michael Shirts',
    author_email='ctk3b@virginia.edu, ctl4f@virginia.edu, edz3fz@virginia.edu, '
                 'michael.shirts@virginia.edu',
    url='https://github.com/shirtsgroup/intermol',
    download_url='https://github.com/shirtsgroup/intermol/tarball/{}'.format(
        intermol.version.short_version),
    packages=find_packages(),
    package_dir={'intermol': 'intermol'},
    include_package_data=True,
    install_requires=reqs,
    license="MIT",
    zip_safe=False,
    keywords='intermol',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
    test_suite='intermol.tests',
    cmdclass={
      'tests': PyTest,
    },
)
