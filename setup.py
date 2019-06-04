import sys

if not sys.version_info[0] == 2:
    print("Invalid version of python, use python 2.")
    sys.exit(1)

import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, Extension

setup(name="reveal", author="Jasper Linthorst", author_email="jasper.linthorst@gmail.com", version="0.2.1",
        url="https://github.com/jasperlinthorst/reveal",
        install_requires=['intervaltree','networkx==2','pysam'],
        packages = ['reveal'],
        scripts = ['ez_setup.py'],
        test_suite = 'nose.collector',
        tests_require = ['nose'],
        ext_modules=[ \
                
                Extension("reveallib", ["reveallib/reveal.c","reveallib/interface.c","divsufsort/divsufsort.c","divsufsort/utils.c","divsufsort/sssort.c","divsufsort/trsort.c"], \
                                       include_dirs=['reveallib','divsufsort'], \
                                       libraries=['pthread'], \
                                       define_macros=[('HAVE_CONFIG_H',1),('__STDC_CONSTANT_MACROS',1),('__STDC_FORMAT_MACROS',1),('__STDC_LIMIT_MACROS',1)], \
                                       undef_macros=['NDEBUG'] ), \
                
                Extension("reveallib64", ["reveallib/reveal.c","reveallib/interface.c","divsufsort/divsufsort.c","divsufsort/utils.c","divsufsort/sssort.c","divsufsort/trsort.c"], \
                                       include_dirs=['reveallib','divsufsort'], \
                                       libraries=['pthread'], \
                                       define_macros = [('SA64',1),('BUILD_DIVSUFSORT64',1),('HAVE_CONFIG_H',1),('__STDC_CONSTANT_MACROS',1),('__STDC_FORMAT_MACROS',1), ('__STDC_LIMIT_MACROS',1)], \
                                       undef_macros=['NDEBUG'] ), \
                
                Extension("probconslib", ["probcons/Probcons.cc"], \
                                       include_dirs=['probcons'], \
                                       define_macros=[('NumInsertStates',2),('VERSION',1.12)], \
                                       undef_macros=['NDEBUG'] ), \
                ],
        entry_points = {
        'console_scripts': [
            'reveal = reveal.reveal:main'
            ]
         }
        )

