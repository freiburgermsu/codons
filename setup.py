# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'Codons',      
  package_dir = {'genes':'codons'},
  packages = find_packages(),
  package_data = {
	'codons':['rosetta_stone/*'],
  },
  version = '0.0.1',
  license = 'MIT',
  description = "Translates and transcribes an arbitrary genetic sequence, and interfaces with BLAST databases to identify the genetic and protein sequences from FASTA formatted files.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/codons',   
  keywords = ['chemistry', 'biology', 'dogma', 'nucleic', 'acids', 'amino', 'translation', 'molecular', 'genetics'],
#  install_requires = ['chemicals', 'pandas']
)