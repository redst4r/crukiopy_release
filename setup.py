from setuptools import setup, find_packages

setup(name='crukiopy_release',
      version=0.1,
      description='Code for the STORMing Cancer Esophagus Atlas manuscript',
      url='http://github.com/redst4r/crukiopy_release/',
      author='Michael Strasser',
      maintainer='Michael Strasser',
      maintainer_email='mstrasse@isbscience.org',
      license='GNU GPL 3',
      keywords='scanpy, scrnaseq',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'sctools @git+https://github.com/redst4r/sctools',  # should pull in all other deps
          'rnaseqtools @git+https://github.com/redst4r/rnaseqtools', 
          'tidyverse @git+https://github.com/redst4r/pytidyverse', 
          'pybustools @git+https://github.com/redst4r/pybustools',
          'sccoda',
          ],
      zip_safe=False
)
