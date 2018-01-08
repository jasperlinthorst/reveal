import sys

if not sys.version_info[0] == 2:
    print("Invalid version of python, use python 2.")
    sys.exit(1)

import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, Extension

setup(name="reveal", version="0.1",
        install_requires=['intervaltree','networkx==2'],
        packages = ['reveal'],
        test_suite = 'nose.collector',
        tests_require = ['nose'],
        ext_modules=[ \
                #Extension("reveallib", ["reveal.c","interface.c","divsufsort/divsufsort.c","divsufsort/utils.c","divsufsort/sssort.c","divsufsort/trsort.c"], \
                #                         libraries=['z','pthread'], \
                #                         #define_macros=[('REVEALDEBUG',1),('_OPENMP',1),('HAVE_CONFIG_H',1),('__STDC_CONSTANT_MACROS',1),('__STDC_FORMAT_MACROS',1),('__STDC_LIMIT_MACROS',1)], \
                #                         define_macros=[('REVEALDEBUG',1),('HAVE_CONFIG_H',1),('__STDC_CONSTANT_MACROS',1),('__STDC_FORMAT_MACROS',1),('__STDC_LIMIT_MACROS',1)], \
                #                         undef_macros=['NDEBUG'] ), \

                Extension("reveallib", ["reveal.c","interface.c","divsufsort/divsufsort.c","divsufsort/utils.c","divsufsort/sssort.c","divsufsort/trsort.c"], \
                                       libraries=['pthread'], \
                                       define_macros=[('HAVE_CONFIG_H',1),('__STDC_CONSTANT_MACROS',1),('__STDC_FORMAT_MACROS',1),('__STDC_LIMIT_MACROS',1)], \
                                       undef_macros=['NDEBUG'] ), \

                Extension("reveallib64", ["reveal.c","interface.c","divsufsort/divsufsort.c","divsufsort/utils.c","divsufsort/sssort.c","divsufsort/trsort.c"], \
                                        libraries=['pthread'], \
                                        define_macros = [('SA64',1),('BUILD_DIVSUFSORT64',1),('HAVE_CONFIG_H',1),('__STDC_CONSTANT_MACROS',1),('__STDC_FORMAT_MACROS',1), ('__STDC_LIMIT_MACROS',1)], \
                                        undef_macros=['NDEBUG'] ), \
                ],
        entry_points = {
        'console_scripts': [
            'reveal = reveal.reveal:main'
            ]
         }
        )

