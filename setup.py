import subprocess
import os
from setuptools import setup

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