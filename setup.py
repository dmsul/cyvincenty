from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    name='Cythonized vincenty',
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize("cyvincenty.pyx"),
)
