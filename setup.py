#
# setup.py
#
# Installation script to get setuptools to install pmutt into
# a Python environment.
#

import sys
import setuptools

# Import the lengthy rich-text README as the package's long
# description:
with open('README.rst', 'r') as fh:
	long_description = fh.read()

setuptools_info = {
	'name': 'pmutt',
	'version': '1.3.0',
	'author': 'Vlachos Research Group',
	'author_email': 'vlachos@udel.edu',
	'description': 'Python Multiscale Thermochemistry Toolbox (pmutt)',
	'long_description': long_description,
	'zip_safe': False,
	'url': 'https://github.com/VlachosGroup/pmutt',
	'packages': setuptools.find_packages(),
	'package_data': {'':['*.xlsx', '*.log', '*OUTCAR']},
	'install_requires': [
		'ASE>=3.16.2',
		'matplotlib>=2.2.3',
		'numpy>=1.15.1',
		'scipy>=1.1.0',
		'pandas>=0.24.2',
		'pymongo>=3.7.2',
		'dnspython>=1.16.0',
		'networkx>=2.3',
		'pygal>=2.4.0',
		'xlrd>=1.2.0',
		'more_itertools>=7.2.0',
	    'PyYAML>=5.1.2'],
	'classifiers': [
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
		"Intended Audience :: Science/Research",
		"Topic :: Scientific/Engineering :: Chemistry",
	    ],
    }

if sys.version_info[0] >= 3:
	#
	# Augment for Python 3 setuptools:
	#
	setuptools_info['long_description_content_type'] = 'text/x-rst'

setuptools.setup(**setuptools_info)
