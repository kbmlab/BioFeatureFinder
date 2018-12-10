from distutils.core import setup
from setuptools import find_packages

setup(
    name='biofeaturefinder',
    version='1.1.4',
    packages=find_packages(),
    url='https://github.com/kbmlab/BioFeatureFinder',
    license='GPLv3',
    author='fciamponi',
    author_email='felipe.ciamponi@gmail.com',
    description='Methods for extracting features from genomic ranges and determining distinguishing marks.',
    install_requires=[
        'glob2 >= 0.6',
        'matplotlib >= 2.2.3',
        'numpy >= 1.15.4',
        'pandas >= 0.23.4',
        'pybedtools >= 0.8.0',
        'pysam >= 0.15.0',
        "rpy2 >= 2.9.4 ; python_version > '3'",
        "rpy2 ; python_version < '3'",
        'scipy >= 1.1.0',
        'seaborn >= 0.9.0',
        'scikit-learn >= 0.20',
    ],
    scripts=[
        'biofeatures/scripts/analyze_features.py',
        'biofeatures/scripts/build_datamatrix.py',
	'biofeatures/scripts/extract_gtf_regions.py',
	'biofeatures/scripts/analyze_gtf_regions.py'
    ]
)
