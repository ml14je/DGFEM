#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joseph Elmes: NERC-Funded PhD Researcher in Applied Mathematics
University of Leeds : Leeds LS2 9JT : ml14je@leeds.ac.uk
Python 3.8: Wed Feb  25 14:03:58 2022
"""

from setuptools import setup

setup(name='DGFEM',
      version= '0.1',
      description= 'A Python translation of DG-FEM code covered in Warburton and Hesthaven, 2008.',
      url= 'http://github.com/ml14je/DGFEM',
      author= 'Joseph Elmes',
      author_email= 'ml14je@leeds.ac.uk',
      license='None',
      install_requires=[
          'wheel', 'numpy', 'scipy', 'pandas', 'matplotlib', 'opencv-python', 'pillow', 'pytube', 'sympy', 'netCDF4', 'plotly'
      ],
      zip_safe=False)
