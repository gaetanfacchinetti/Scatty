import subprocess
import os

from setuptools import setup, Extension
from Cython.Build import cythonize

"""
def main():

    # go to where the project needs to be compiled
    
    os.chdir("./src/scatty/src/")

    try:

        if os.path.exists("libcscatty.so"):
            os.remove("libcscatty.so")

        subprocess.run(['make'], check=True)
        print("make executed successfully.")

        # Run the 'make' command with the desired target
        subprocess.run(['make', 'clean'], check=True)
        print("make clean executed successfully.")


    except subprocess.CalledProcessError as e:
        print("Error executing Makefile:", e)

    # go back to the main project folder
    os.chdir("../../../")
    
    # cytonize the cython code
    cythonize("./src/scatty/wrapper.pyx", force = True)

    # define the extension for the libraries
    extension = Extension("scatty.wrapper", 
                        sources = ["./src/scatty/wrapper.c"], 
                        libraries = ["cscatty", "m"], 
                        library_dirs=["./src/scatty/src/"],)
    
    setup(ext_modules=[extension], 
          package_dir={'' : 'src'},  
          packages=['scatty'], 
          package_data= {'scatty' : ['data/*', 'src/*.h', 'src/*.so']},
          include_package_data=True)


if __name__ == "__main__":
    main()
"""

from distutils.command.install import install as _install


class install(_install):
    def run(self):
        os.chdir("./src/scatty/src/")
        subprocess.run(['make', 'clean'], check=True)
        subprocess.run(['make', 'cleanlib'], check=True)
        subprocess.run(['make'], check=True)
        subprocess.run(['make', 'clean'], check=True)
        os.chdir("../../../")
        _install.run(self)


setup(
    name='scatty',
    version='0.0.1',
    package_dir={'' : 'src'},
    packages=['scatty'],
    package_data= {'scatty' : ['data/*', 'src/*.h', 'src/*.so']},
    include_package_data=True,
    cmdclass={'install': install},
)