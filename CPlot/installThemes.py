# install les themes ttk
import KCore.installPath
import os, sys
force = False
if len(sys.argv) > 1:
    if sys.argv[1] == '-f': force = True
    elif sys.argv[1] == '-force': force = True

p = KCore.installPath.installPath+'/CPlot'
a = os.access(p+'/themes.tar', os.R_OK)
b = os.access(p+'/themes', os.R_OK)
if force or (a and not b):
    h = os.getcwd()
    os.chdir(p)
    os.system('tar xvf themes.tar')
    os.chdir(h)
