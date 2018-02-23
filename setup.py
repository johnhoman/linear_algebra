from distutils.core import setup, Extension

extension_mod = Extension("c_gaussian", ["gaussian.c"], undef_macros=['NDEBUG'])
setup(name="c_gaussian", ext_modules=[extension_mod])