from distutils.core import setup
from setuptools import find_packages

setup(
    name='biofeaturefinder',
    version='1.1.3',
    packages=find_packages(),
    url='https://github.com/kbmlab/BioFeatureFinder',
    license='GPLv3',
    author='fciamponi',
    author_email='felipe.ciamponi@gmail.com',
    description='Methods for extracting features from genomic ranges and determining distinguishing marks.',
    install_requires=[
#        'argparse',
	'glob2',
	'matplotlib',
	'numpy',
	'pandas',
	'pybedtools',
	'pysam',
	'rpy2',
	'scipy',
	'seaborn',
	'scikit-learn',
#	'system'
    ],
    scripts=[
        'biofeatures/scripts/analyze_features.py',
        'biofeatures/scripts/build_datamatrix.py',
	'biofeatures/scripts/extract_gtf_regions.py',
	'biofeatures/scripts/analyze_gtf_regions.py'
    ]
)
