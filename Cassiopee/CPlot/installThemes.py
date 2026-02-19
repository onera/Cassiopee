# install les themes ttk
import KCore.installPath as K
import os
import sys

force = False
if len(sys.argv) > 1:
    if sys.argv[1] in ['-f', '-force']: force = True

p = os.path.join(K.installPath, 'CPlot')
a = os.access(os.path.join(p, 'themes.tar'), os.R_OK)
b = os.access(os.path.join(p, 'themes'), os.R_OK)
if force or (a and not b):
    h = os.getcwd()
    os.chdir(p)
    os.system('tar xvf themes.tar')
    os.chdir(h)
