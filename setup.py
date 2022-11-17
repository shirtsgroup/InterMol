"""InterMol: A conversion tool for molecular dynamics simulations.
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = None

setup(
    name='intermol',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=short_description,
    long_description=long_description,
    author='Christoph Klein, Christopher Lee, Ellen Zhong, and Michael Shirts',
    author_email='ctk3b@virginia.edu, ctl4f@virginia.edu, edz3fz@virginia.edu, '
                 'michael.shirts@virginia.edu',
    url='https://github.com/shirtsgroup/intermol',
    packages=find_packages(),
    package_dir={'intermol': 'intermol'},
    package_data={'tests': ['*.py',
                            '*.md',
                            'desmond/*.cfg',
                            'desmond/*/*.cms',
                            'gromacs/*.mdp',
                            'gromacs/*/*/*.gro',
                            'gromacs/*/*/*.top',
                            'gromacs/*/*/*.itp',
                            'lammps/*/*.lmp',
                            'lammps/*/*.input',
                            'amber/*.in',
                            'amber/*/*',
                            ]},
    include_package_data=True,
    license="MIT",
    zip_safe=False,
    keywords='intermol',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
    entry_points={
        'console_scripts': [
            'intermol-convert=intermol.convert:main',
        ],
    },
)

