#=============================================================================
# Apps requires:
# ELSAPROD variable defined in environment
# CASSIOPEE
#=============================================================================
import os
from setuptools import setup
from importlib.util import spec_from_file_location, module_from_spec

def loadModuleFromPath(modname):
    # Load a Python file by filesystem path (PEP-517 isolated build requirement)
    helper = os.path.join(os.path.dirname(__file__), modname + ".py")
    spec = spec_from_file_location(modname, helper)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# Compiler settings must be set in installBase.py / installBaseUser.py
Dist = loadModuleFromPath('../KCore/Dist')

prod = os.getenv("ELSAPROD") or "xx"

# setup ======================================================================
setup(
    name="Apps",
    version="4.1",
    description="Application modules",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Apps', 'Apps.Chimera', 'Apps.Fast', 'Apps.Mesh', 'Apps.Coda'],
    package_dir={"":"."}
)

# Check PYTHONPATH ===========================================================
installPath = loadModuleFromPath('../KCore/installPath')
installPathDict = {
    "installPath": installPath.installPath,
    "libPath": installPath.libPath,
    "includePath": installPath.includePath
}
Dist.checkPythonPath(installPathDict)
Dist.checkLdLibraryPath(installPathDict)
