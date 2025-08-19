#from distutils.core import setup # Seems like obsolete
from setuptools import setup
from Cython.Build import cythonize
setup(ext_modules = cythonize('gen_ngrams_cy.py'))
