from distutils.core import setup
from distutils.sysconfig import get_python_lib
import glob
import os
import sys

setup(
    name = "Pypdm",
    packages     = ['Pypdm'],
    data_files = [('', ["Pypdm.so"])],
    author = 'E. Quemerais',
    description = 'Toolkit for parallel distributed computational geometry',
    license = 'LGPL',
    )

