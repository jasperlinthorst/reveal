import sys

if not sys.version_info[0] == 2:
    print("Invalid version of python, use python 2.")
    sys.exit(1)

import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, Extension

setup(name="reveal", version="0.1",
        install_requires=['intervaltree','networkx'],
        packages = ['reveal'],
        ext_modules=[ \
                Extension("reveallib", ["reveal.c","interface.c"], libraries=['z','divsufsort','pthread'], define_macros = [], undef_macros=['NDEBUG'] ), \
                Extension("reveallib64", ["reveal.c","interface.c"], libraries=['z','divsufsort64','pthread'], define_macros = [('SA64',1)], undef_macros=['NDEBUG'] ), \
                ],
        entry_points = {
        'console_scripts': [
            'reveal = reveal.reveal:main'
            ]
         }
        )

