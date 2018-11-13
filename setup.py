"""
Qcat setup script

Copyright (c) 2018 by Oxford Nanopore Technologies Ltd.
"""
from setuptools import setup

from qcat import __version__

setup(
    name='qcat',
    version=__version__,
    url='https://github.com/nanoporetech/qcat/',
    license='Mozilla Public License Version 2.0',
    description="Qcat is Python command-line tool for demultiplexing Oxford Nanopore reads from FASTQ files",
    author='Oxford Nanopore Technologies Ltd.',
    author_email='philipp.rescheneder@nanoporetech.com',
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
        'pytest-runner',
        'pytest',
        'pytest-cov',
    ],
    packages=[
        'qcat'
    ],
    package_data={
        'qcat': [
            'resources/kits/*',
        ]
    },

    entry_points={"console_scripts": ['qcat = qcat.cli:main',
                                      'qcat-eval = qcat.eval:main',
                                      'qcat-roc = qcat.eval_roc:main',
                                      'qcat-eval-truth = qcat.eval_full:main']}
)