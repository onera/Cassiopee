# Installation de libkcore pour etre accessible par les autres modules
# Si libkcore.a existe, on la recopie
# Sinon, on cherche kcore.so ou kcore.pyd, on le recopie en libkcore.so ou dll
# Installe aussi les .py qui ne sont pas installes par install
# Installe aussi les fichiers d'environnement pour l'execution
import os, shutil
import Dist

# Symlinks eventuel
Dist.symLinks()

system = Dist.getSystem()[0]

if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

try: import KCore.installPath as K
except: import installPath as K
libPath = K.libPath
installPathLocal = K.installPath

# La librarie statique existe?
a = os.access(installPathLocal+"/KCore/libkcore.a", os.F_OK)
if a:
    shutil.copy(installPathLocal+"/KCore/libkcore.a", libPath+"/libkcore.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/KCore/kcore"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copy(installPathLocal+"/KCore/kcore"+__EXTMODULE__,
                    libPath+"/libkcore"+__EXTSHARED__) 
    else:
        print "Error: kcore%s can not be found in %s."%(__EXTMODULE__,installPathLocal)

# Copie aussi les .py
shutil.copy("config.py", installPathLocal+"/KCore/config.py")
shutil.copy("Dist.py", installPathLocal+"/KCore/Dist.py")
shutil.copy("installPath.py", installPathLocal+"/KCore/installPath.py")
shutil.copy("installBase.py", installPathLocal+"/KCore/installBase.py")

# Ecrit les infos d'install
import Dist
Dist.writeBuildInfo()
shutil.copy("buildInfo.py", installPathLocal+"/KCore/buildInfo.py")

# Ecrit les fichiers d'environnement
Dist.writeEnvs()

# Installe la licence
f = open('KCore/installKey.py'); a = f.read(); f.close()
exec a
