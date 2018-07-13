# install les themes ttk
import KCore.installPath
import os
p = KCore.installPath.installPath+'/CPlot'
a = os.access(p+'/themes.tar', os.R_OK)
b = os.access(p+'/themes', os.R_OK)
if a and not b:
    h = os.getcwd()
    os.chdir(p)
    os.system('tar xvf themes.tar')
    os.chdir(h)
