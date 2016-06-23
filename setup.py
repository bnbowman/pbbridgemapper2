from setuptools import find_packages, setup

import subprocess
import sys

version = __import__('pbbridgemapper').get_version()

returnCode = subprocess.call('blasr', shell=True)
if returnCode != 0:
    sys.stderr.write("Unable to install -- must install blasr first.\n")
    sys.exit(1)

setup(
    name='pbbridgemapper',
    version=version,
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    license='LICENSE.txt',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'pbbridgemapper = pbbridgemapper.main:main']
    },
    zip_safe=False,
    install_requires=['pbcore >= 0.6.0']
)
