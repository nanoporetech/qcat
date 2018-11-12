"""
Qcat setup script

Copyright (c) 2018 by Oxford Nanopore Technologies Ltd.
"""
from setuptools import setup

from qcat import __version__

setup(
    name='qcat',
    version=__version__,
    url='',
    author='Oxford Nanopore Technologies Ltd.',
    author_email='philipp.rescheneder@nanoporetech.com',
    description='Demultiplexing Nanopore read FASTQ files',
    setup_requires=['pytest-runner'],
    install_requires=[
        'biopython',
        'parasail',
        'six',
        'pyyaml'
    ],
    extras_require={
        'DOCS': ['sphinx>=1.3.1'],
    },
    tests_require=[
        'pytest',
        'pytest-cov',
    ],
    packages=[
        'qcat'
    ],
    package_data={
        'qcat': [
            'resources/*',
            'resources/kits/*',
            'resources/models/*',
        ]
    },

    entry_points={"console_scripts": ['qcat = qcat.cli:main',
                                      'qcat-eval = qcat.eval:main',
                                      'qcat-roc = qcat.eval_roc:main',
                                      'qcat-eval-truth = qcat.eval_full:main']}
)