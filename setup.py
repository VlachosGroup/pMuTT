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
	'version': '1.4.1',
	'author': 'Vlachos Research Group',
	'author_email': 'vlachos@udel.edu',
	'description': 'Python Multiscale Thermochemistry Toolbox (pmutt)',
	'long_description': long_description,
	'zip_safe': False,
	'url': 'https://github.com/VlachosGroup/pmutt',
	'packages': setuptools.find_packages(),
	'package_data': {'':['*.xlsx', '*.log', '*OUTCAR']},
	'install_requires': [
		'ASE>=3.22.1',
		'matplotlib>=3.6.2',
		'numpy>=1.21.5',
		'scipy>=1.9.3',
		'pandas>=1.4.4',
		'pymongo>=4.2.0',
		'dnspython>=2.2.1',
		'networkx>=2.8.4',
		'pygal>=3.0.0',
		'xlrd>=2.0.1',
		'more_itertools>=8.14.0',
	    'PyYAML>=6.0.0',
		'python>=3.3'
		],
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
