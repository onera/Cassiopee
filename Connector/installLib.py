# Installation de libconnector pour etre accessible par les autres modules
# Si libconnector.a existe, on la recopie
# Sinon, on cherche connector.so ou connector.pyd, on le recopie en libconnector.so ou dll
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
a = os.access(installPathLocal+"/Connector/libconnector.a", os.F_OK)
if a:
    shutil.copy(installPathLocal+"/Connector/libconnector.a", libPath+"/libconnector.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/Connector/connector"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copy(installPathLocal+"/Connector/connector"+__EXTMODULE__,
                    libPath+"/libconnector"+__EXTSHARED__) 
    else:
        print("Error: connector"+__EXTMODULE__+" can not be found.")
