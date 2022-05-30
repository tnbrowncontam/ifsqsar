"""
ifsqsar/setup.py
developed by Trevor N. Brown
setup file for ifsqsar
"""

import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()

setuptools.setup(
    name='ifsqsar',
    version='1.0.1',
    author='Trevor N. Brown',
    author_email='trevor.n.brown@gmail.com',
    description='A package for applying IFS QSARs',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/tnbrowncontam/ifsapp',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: OS Independent',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    python_requires='>=3.4',
)

