# This file must replace installPath.py when producing AppImage
import os
HERE = os.environ['APPDIR']
installPath = HERE+'/Dist/lib/python2.7/site-packages'
libPath = HERE+'/Dist/lib'
includePath = ''
