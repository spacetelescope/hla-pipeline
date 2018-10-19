#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', 'astropy', 'scipy', 'matplotlib', \
                'photutils', 'pyyaml', 'astroquery',\
                'scipy', 'sphinx','sphinx_rtd_theme', \
                'stsci_rtd_theme','stsci.tools',\
                'stwcs','setuptools']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'requests_mock', 'ci_watson']

setup(
    author="Warren J. Hack",
    author_email='hack@stsci.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Code for implementing HLA-type processing in the HST Pipeline",
    entry_points={
        'console_scripts': [
            'hlapipeline=hlapipeline.cli:main',
        ],
    },
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='hlapipeline',
    name='hlapipeline',
    packages=find_packages(),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/spacetelescope/hlapipeline',
    version='0.1.0',
    zip_safe=False,
)
