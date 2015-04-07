from distutils.core import setup, Extension
setup(name="GSA", version="1.0",
      ext_modules=[ \
			Extension("GSA_64", ["GSA.c"], libraries=['z','divsufsort64'], undef_macros=['NDEBUG'], define_macros = [('LARGEIDX', '1')], extra_compile_args=['-O3']), \
			Extension("GSA", ["GSA.c"], libraries=['z','divsufsort'], undef_macros=['NDEBUG'], extra_compile_args=['-O3']), \
			])

