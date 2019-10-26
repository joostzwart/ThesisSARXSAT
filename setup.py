# Created by Joost Zwart as part of the master thesis project at the TU Delft
#
# Still under heavy development. Use at your own risk.
import numpy
from setuptools import setup
from distutils.extension import Extension

setup(  name='SarxIdent',
        version='0.1',
        description='Identification tool for switched ARX systems. UNDER CONSTRUCTION, use at your own risk',
        author='Joost Zwart',
        install_requires=['z3-solver','cplex','PyQt5','numpy','matplotlib'],
        ext_modules = [Extension("block_operations", ["build/pythonfiles/block_operations.c"],include_dirs=[numpy.get_include()]),
                       Extension("checkPoints", ["build/pythonfiles/checkPoints.c"], include_dirs=[numpy.get_include()]),
                       Extension("DataGenerator", ["build/pythonfiles/DataGenerator.c"],include_dirs=[numpy.get_include()]),
                       Extension("dataset", ["build/pythonfiles/dataset.c"],include_dirs=[numpy.get_include()]),
                       Extension("export",["build/pythonfiles/export.c"],include_dirs=[numpy.get_include()]),
                       Extension("identify",["build/pythonfiles/identify.c"],include_dirs=[numpy.get_include()]),
                       Extension("L0",["build/pythonfiles/L0.c"],include_dirs=[numpy.get_include()]),
                       Extension("LinProg",["build/pythonfiles/LinProg.c"],include_dirs=[numpy.get_include()]),
                       Extension("model",["build/pythonfiles/model.c"],include_dirs=[numpy.get_include()]),
                       Extension("parameters",["build/pythonfiles/parameters.c"],include_dirs=[numpy.get_include()]),
                       Extension("partition",["build/pythonfiles/partition.c"], include_dirs=[numpy.get_include()]),
                       Extension("QuadProg",["build/pythonfiles/QuadProg.c"], include_dirs=[numpy.get_include()]),
                       Extension("models",["build/pythonfiles/models.c"],include_dirs=[numpy.get_include()]),
                       Extension("ReaderData",["build/pythonfiles/ReaderData.c"],include_dirs=[numpy.get_include()]),
                       Extension("SAT",["build/pythonfiles/SAT.c"],include_dirs=[numpy.get_include()]),
                       Extension("SimplifyCertificate",["build/pythonfiles/SimplifyCertificate.c"],include_dirs=[numpy.get_include()]),
                       Extension("Sparsification",["build/pythonfiles/Sparsification.c"],include_dirs=[numpy.get_include()]),
                       Extension("storage",["build/pythonfiles/storage.c"],include_dirs=[numpy.get_include()]),
                       Extension("TheorySolver",["build/pythonfiles/TheorySolver.c"],include_dirs=[numpy.get_include()]),
                       ]
    )
