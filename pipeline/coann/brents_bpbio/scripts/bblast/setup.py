from setuptools import setup, find_packages
import sys, os


version = '0.1'

setup(name='bblast',
      version=version,
      description="""a dropin replacement for blast that will save params, run
      formatdb, and not re-run if input files and prams are unchanged""",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bioinformatics blast',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests'] +
                             ["."]),
      include_package_data=True,
      test_suite='nose.collector',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=["bblast/bblast.py"],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
