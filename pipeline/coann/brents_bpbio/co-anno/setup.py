from setuptools import setup, find_packages

version = '0.2'

setup(name='coanno',
      version=version,
      description="""Given a pair of gff files, and a pair of fasta files, use
      one to annotate the other for missed exons""",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='bioinformatics',
      author='brentp',
      author_email='bpederse@gmail.com',
      url='',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
          'numpy', 'bblast', 'simplejson'
      ],
      scripts=["coanno/coannotate.py", "coanno/mask_features.py"],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
