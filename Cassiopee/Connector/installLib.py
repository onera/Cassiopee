# Installation de libconnector pour etre accessible par les autres modules
# Si libconnector.a existe, on la recopie
# Sinon, on cherche connector.so ou connector.pyd, on le recopie en libconnector.so ou dll
import os, shutil
import KCore.Dist as Dist
system = Dist.getSystem()[0]

if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

import KCore.installPath as K
libPath = K.libPath
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
installPathLocal = 'build/'+prod

# La librarie statique existe?
a = os.access(installPathLocal+"/libconnector.a", os.F_OK)
if a:
    shutil.copyfile(installPathLocal+"/libconnector.a", libPath+"/libconnector.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/connector"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copyfile(installPathLocal+"/connector"+__EXTMODULE__,
                        libPath+"/libconnector"+__EXTSHARED__)
    else:
        print("Error: connector"+__EXTMODULE__+" can not be found.")
