# coding: utf-8
__author__ = "RCAS development team"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Dilmurat Yusuf", "Bora Uyar",
             "Ricardo Wurmus", "Altuna Akalin"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Dilmurat Yusuf"
__email__ = "dilmurat.yusuf@gmail.com"

long_description = """
RCAS: RNA Centric Annotation System provides intuitive reports
and publication ready graphics.

https://github.com/BIMSBbioinfo/RCAS.git
"""

classifiers = [
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3",
    "Operating System :: Unix",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics"]

import sys

try:
    from setuptools import setup
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.", file=sys.stderr)

setup(name='RCAS',
      version=__version__,
      description='RNA Centric Annotation System',
      long_description=long_description,
      author=__author__,
      classifiers=classifiers,
      author_email=__email__,
      maintainer=___maintainer__,
      maintainer_email=__email__,
      url='https://github.com/BIMSBbioinfo/RCAS.git',
      packages=['RCAS'],
      package_data={'RCAS':
                    ['data/*gmt',
                     'data/*meme',
                     'data/*css',
                     'data/*html',
                     'data/img/*',
                     'data/snakefiles/*'],
                    },
      license=__license__,
      keywords=['bioinformatics', 'microbiome', 'microbiology', 'RCAS'],
      platforms=['Linux'],
      install_requires=['numpy >= 1.9.0',
                        'scipy >= 0.14.0',
                        'cogent == 1.5.3',
                        'natsort < 4.0.0',
                        'matplotlib >= 1.1.0, != 1.4.2',
                        'pynast == 1.2.2', 'qcli >= 0.1.1, < 0.2.0', 'gdata',
                        'biom-format >= 2.1.4, < 2.2.0',
                        'emperor >= 0.9.51, < 1.0.0',
                        'scikit-bio >= 0.2.3, < 0.3.0',
                        'burrito-fillings >= 0.1.1, < 0.2.0',
                        'pandas >= 0.13.1', 'burrito >= 0.9.1, < 1.0.0',
                        'RCAS-default-reference >= 0.1.2, < 0.2.0'],
      extras_require={'all': ['ipython[notebook] >= 3.1.0, < 4.0.0',
                              'sphinx >= 0.3']}
      )
