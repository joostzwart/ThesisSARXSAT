# Created by Joost Zwart as part of the master thesis project at the TU Delft
#
# Still under heavy development. Use at your own risk.
import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
extensions = [

    Extension("*", ["pythonfiles/*.pyx"],
        include_dirs=["/pythonfiles",numpy.get_include()],)
    ]
setup(
    name="SSHIP",
    ext_modules=cythonize(extensions,build_dir="build"),
)

