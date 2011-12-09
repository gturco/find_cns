from setuptools import setup, find_packages
from distutils.extension import Extension
#from Cython.Distutils import build_ext

version = '0.0'

setup(name='biostuff',
      version=version,
      description="",
      long_description="""\
""",
      ext_modules=[ Extension("biostuff/cblastline",
                      sources=["biostuff/cblastline.pyx"],)],
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bio',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='BSD',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      test_suite='nose.collector',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points= {
          # -*- Entry points: -*-
          'console_scripts': ['biostuff = biostuff:main']
          }
      )
