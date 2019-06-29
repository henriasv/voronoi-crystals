import setuptools
from setuptools import setup

setup(name='voronoi-crystals',
      version='0.1',
      description='Tool to create polycrystals by voronoi tesselation',
      url='http://github.com/henriasv/voronoi-crystals',
      author='Henrik Andersen Sveinsson',
      author_email='henrik.sveinsson@me.com',
      license='GNU GPL v3.0',
      packages=setuptools.find_packages(),
      zip_safe=False, 
      include_package_data=True)