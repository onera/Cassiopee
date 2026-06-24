# install les themes ttk
import KCore.installPath as K
import os
import sys
import subprocess

force = False
if len(sys.argv) > 1:
    if sys.argv[1] in ['-f', '-force']: force = True

p = os.path.join(K.installPath, 'CPlot')
a = os.access(os.path.join(p, 'themes.tar'), os.R_OK)
b = os.access(os.path.join(p, 'themes'), os.R_OK)
if force or (a and not b):
    subprocess.run(["tar", "xvf", "themes.tar"], cwd=p, check=False)
