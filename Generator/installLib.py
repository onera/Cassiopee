# Installation de libgenerator pour etre accessible par les autres modules
# Si libgenerator.a existe, on la recopie
# Sinon, on cherche generator.so ou generator.pyd, on le recopie en 
# libgenerator.so ou dll
import os, shutil
import platform
system = platform.uname()[0]

if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

import KCore.installPath as K
libPath = K.libPath
installPathLocal = K.installPath

# La librarie statique existe?
a = os.access(installPathLocal+"/Generator/libgenerator.a", os.F_OK)
if a:
    shutil.copy(installPathLocal+"/Generator/libgenerator.a", libPath+"/libgenerator.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/Generator/generator"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copy(installPathLocal+"/Generator/generator"+__EXTMODULE__,
                    libPath+"/libgenerator"+__EXTSHARED__) 
    else:
        print("Error: generator"+__EXTMODULE__+" can not be found.")
