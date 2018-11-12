from distutils.core import setup
from setuptools import find_packages

setup(
    name='featurefinder',
    version='1.1.1',
    packages=find_packages(),
    url='https://github.com/kbmlab/BioFeatureFinder',
    license='',
    author='fciamponi',
    author_email='felipe.ciamponi@gmail.com',
    description='Methods for extracting features from genomic ranges and determining distinguishing marks.',
    install_requires=[
        'pandas',
        'scikit-learn',
        'pybedtools',
        'pysam',
        'matplotlib',
        'scipy',
        'rpy2'
    ],
    scripts=[
        'biofeatures/scripts/analyze_features.py',
        'biofeatures/scripts/build_datamatrix.py',
	'biofeatures/scripts/extract_gtf_regions.py',
	'biofeatures/scripts/analyze_gtf_regions.py'
    ]
)
