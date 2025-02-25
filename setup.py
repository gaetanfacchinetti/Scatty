from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

# Ensure absolute paths
include_gsl_dir = "/usr/local/include"
lib_gsl_dir = "/usr/local/lib"

ext_modules = [
    Extension(
        "scatty.libscatty",
        sources=["src_c/cross_section.c", "src_c/integration.c"],
        include_dirs=[numpy.get_include(), include_gsl_dir, "src_c"],
        library_dirs=[lib_gsl_dir],
        libraries=["gsl", "gslcblas", "m"],
        language="c",
        extra_compile_args=["-std=c11"],
    )
]

setup(
    name="scatty",
    version="0.1.0",
    packages=["scatty"],
    ext_modules=cythonize(ext_modules, compiler_directives={"language_level": "3"}),
)