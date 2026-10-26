# Installation de libxcore pour etre accessible par les autres modules
# Si libxcore.a existe, on la recopie
# Sinon, on cherche xcore.so ou xcore.pyd, on le recopie en
# libxcore.so ou dll
import os
import shutil
import KCore.Dist as Dist

system = Dist.getSystem()[0]
if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

libPath = Dist.getLibPath()
prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")
installPathLocal = 'build/'+prod

# La librarie statique existe?
a = os.access(installPathLocal+"/libxcore.a", os.F_OK)
if a:
    shutil.copyfile(installPathLocal+"/libxcore.a", libPath+"/libxcore.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/xcore"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copyfile(installPathLocal+"/xcore"+__EXTMODULE__,
                        libPath+"/libxcore"+__EXTSHARED__)
    else:
        print("Error: xcore"+__EXTMODULE__+" can not be found.")
