#from distutils.core import setup, Extension
from setuptools import setup, Extension

setup(name="reveal", version="0.1",
	install_requires=['seqalign==0.1',],
	dependency_links = ['https://github.com/mhulsman/seqal/archive/seqal_python_v0.1.tar.gz#egg=seqal-0.1'],
	ext_modules=[ \
		Extension("GSA_64", ["GSA.c"], libraries=['z','divsufsort64'], undef_macros=['NDEBUG'], define_macros = [('LARGEIDX', '1')], extra_compile_args=['-O3']), \
		Extension("GSA", ["GSA.c"], libraries=['z','divsufsort'], undef_macros=['NDEBUG'], extra_compile_args=['-O3']), \
	]
)
