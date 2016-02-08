# coding: utf-8
__author__ = "RCAS development team"
__copyright__ = "Copyright (c) 2016--, %s" % __author__
__credits__ = ["Dilmurat Yusuf", "Bora Uyar",
             "Ricardo Wurmus", "Altuna Akalin"]
__license__ = "MIT"
__version__ = "0.1.0"
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
    "Programming Language :: Python :: 2.7",
    "Operating System :: Unix",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics"]

from setuptools import setup

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
      license=__license__,
      keywords=['Bioinformatics', 'Clip-Seq', 'Peaks', 'Annotation'],
      platforms=['Linux'],
      entry_points={"console_scripts": ["RCAS = RCAS.RCAS:main"]}
      package_data={'RCAS':
                    ['data/*gmt',
                     'data/*meme',
                     'data/*css',
                     'data/*html',
                     'data/img/*',
                     'data/snakefiles/*',
                     'libexec/*R',
                     'libexec/*py',
                     'libexec/rcas.Rmd',
                     'libexec/generate_report.sh'],
                    },
      install_requires=[
                      'snakemake',
                      'bedtools < 2.5.0',
                      'MEME-chip',
                      'pandoc >= 1.12.3',
                      'rmarkdown',
                      'rtracklayer',
                      'data.table',
                      'biomaRt',
                      'org.Hs.eg.db',
                      'org.Ce.eg.db',
                      'org.Mm.eg.db',
                      'org.Dm.eg.db',
                      'topGO',
                      'DT',
                      'plotly',
                      'dplyr',
                      'genomation',
                      'GenomicFeatures'
                      ]
      )
