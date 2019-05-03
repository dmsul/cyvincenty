from setuptools import setup, Extension
from Cython.Distutils import build_ext

import numpy as np

EXTENSIONS = [
    Extension('cyvincenty.core',
              ['cyvincenty/core.pyx'],
              include_dirs=[np.get_include()],
              extra_compile_args=['/openmp', '/O2', '/fp:fast'],
              ),
]

setup(
    name='cyvincenty',
    version='0.1',
    packages=['cyvincenty'],
    cmdclass={'build_ext': build_ext},
    ext_modules=EXTENSIONS,
)
