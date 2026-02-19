#=============================================================================
# Apps requires:
# ELSAPROD variable defined in environment
# CASSIOPEE
#=============================================================================
import os
from setuptools import setup
import KCore.Dist as Dist

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

