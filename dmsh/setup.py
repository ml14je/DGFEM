#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joseph Elmes: NERC-Funded PhD Researcher in Applied Mathematics
University of Leeds : Leeds LS2 9JT : ml14je@leeds.ac.uk

Python 3.7
Created on Thu Jan 14 19:42:02 2021
"""
from setuptools import setup

setup(name='dmsh v2 package',
      version= '0.1',
      description= 'Originally written by Nico Schlömer (nico.schloemer@gmail.com), the following code produces *and refines* high quality meshes using the delauney method',
      url= 'http://github.com/ml14je/dmshv2',
      author= 'Joseph Elmes',
      author_email= 'ml14je@leeds.ac.uk',
      license='None',
      install_requires=[
          'meshplex', 'numpy', 'scipy'
      ],
      zip_safe=False)