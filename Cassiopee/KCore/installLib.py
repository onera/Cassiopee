# Installation de libkcore pour etre accessible par les autres modules
# Si libkcore.a existe, on la recopie
# Sinon, on cherche kcore.so ou kcore.pyd, on le recopie en libkcore.so ou dll
# Installe aussi les .py qui ne sont pas installes par install
# Installe aussi les fichiers d'environnement pour l'execution
import os
import sys
import shutil
import Dist

system = Dist.getSystem()[0]
if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

# Check PYTHONPATH
prefix = Dist.parseInstallPathPrefix(*sys.argv)
installPath = os.path.join(prefix, Dist.getInstallPathSuffix())
libPath = os.path.dirname(os.path.dirname(installPath))
Dist.checkPythonPath(installPath)
Dist.checkLdLibraryPath(libPath)

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")
installPathLocal = 'build/'+prod

# La librarie statique existe?
a = os.access(installPathLocal+"/libkcore.a", os.F_OK)
if a:
    shutil.copyfile(installPathLocal+"/libkcore.a", libPath+"/libkcore.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/kcore"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copyfile(installPathLocal+"/kcore"+__EXTMODULE__,
                        libPath+"/libkcore"+__EXTSHARED__)
    else:
        print("Error: kcore%s can not be found in %s."%(__EXTMODULE__,installPathLocal))

# Copie aussi les .py
installPath = os.path.join(installPath, "KCore")
shutil.copyfile("Dist.py", installPath+"/Dist.py")
shutil.copyfile("installBase.py", installPath+"/installBase.py")
if os.access("installBaseUser.py", os.R_OK):
    shutil.copyfile("installBaseUser.py", installPath+"/installBaseUser.py")
shutil.copyfile("test/notify.py", installPath+"/notify.py")

# Ecrit les infos d'install
Dist.writeBuildInfo()
shutil.copyfile("buildInfo.py", installPath+"/buildInfo.py")

# Ecrit les fichiers d'environnement
Dist.writeEnvs(installPath, libPath)
