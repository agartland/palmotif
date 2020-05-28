from setuptools import setup, find_packages
PACKAGES = find_packages()

"""Read the contents of your README file"""
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

opts = dict(name='seqlogo',
            maintainer='Andrew Fiore-Gartland',
            maintainer_email='agartlan@fredhutch.org',
            description='Python package for computing and plotting sequence logos, with output as SVG or matplotlib',
            long_description=(""""""),
            long_description_content_type='text/markdown',
            url='https://github.com/agartland/seqlogo',
            license='MIT',
            author='Andrew Fiore-Gartland',
            author_email='agartlan@fredhutch.org',
            version='0.1',
            packages=PACKAGES
           )

requires = ['matplotlib',
            'numpy>=1.16',
            'pandas>=0.24.2',
            'parasail>=1.1.17',
            'svgwrite']
setup(**opts, install_requires=requires)