#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='bigwig-plotter',
    version='0.1',
    packages=['src'],
    url='',
    license='',
    include_package_data=True,
    author='brianyee',
    author_email='',
    description='plots bigwig density signal',
    package_dir={
        'bigwig-plotter': 'src',
    },
    entry_points = {
        'console_scripts': [
            'bigwig-plotter = src.plotter:main',
        ]
    }
)
