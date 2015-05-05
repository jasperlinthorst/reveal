import sys

if not sys.version_info[0] == 2:
    print("Invalid version of python, use python 2.")
    sys.exit(1)

import ez_setup
ez_setup.use_setuptools()

from setuptools import setup, Extension

setup(name="reveal", version="0.1",
	install_requires=['seqal==0.1',],
	dependency_links = ['https://github.com/mhulsman/seqal/archive/seqal_python_v0.1.tar.gz#egg=seqal-0.1'],
	scripts = ['reveal'],
	ext_modules=[ \
		Extension("GSA_64", ["GSA.c"], libraries=['z','divsufsort64'], undef_macros=['NDEBUG'], define_macros = [('LARGEIDX', '1')], extra_compile_args=['-O3']), \
		Extension("GSA", ["GSA.c"], libraries=['z','divsufsort'], undef_macros=['NDEBUG'], extra_compile_args=['-O3']), \
	]
)
