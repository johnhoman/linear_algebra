from distutils.core import setup, Extension

extension_mod = Extension("c_gaussian", ["row_echelon_transformation.c"], undef_macros=['NDEBUG'])
setup(name="c_gaussian", ext_modules=[extension_mod])
