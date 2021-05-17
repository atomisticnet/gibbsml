from setuptools import setup, find_packages
from codecs import open
import os

__author__ = "Alexander Urban"
__email__ = "aurban@atomistic.net"

here = os.path.abspath(os.path.dirname(__file__))
package_name = 'gibbsml'
package_description = 'Prediction of oxide formation free energies'

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as fp:
    long_description = fp.read()

# Get version number from the VERSION file
with open(os.path.join(here, package_name, 'VERSION')) as fp:
    version = fp.read().strip()

setup(
    name=package_name,
    version=version,
    description=package_description,
    long_description=long_description,
    url='ellingham.energy-materials.org',
    author=__author__,
    author_email=__email__,
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    keywords=['materials science', 'phase diagrams', 'machine learning'],
    packages=find_packages(exclude=['tests']),
    install_requires=['ase>=3.17.0', 'numpy>=1.14.3', 'scipy>=1.1.0',
                      'pandas>=0.24.0', 'pymatgen>=2020.7.3',
                      'scikit-learn>=0.19.1', 'catlearn>=0.6.2',
                      'plotly>=4.2.1', 'mendeleev', 'ase'],
    # package_data={
    #     'sample': [''],
    # },
    entry_points={
        'console_scripts': [
            'dribble=dribble.scripts.dribble:main',
        ],
    },
    include_package_data=True,
    # scripts=glob.glob(os.path.join("scripts", "*.py"))
)
