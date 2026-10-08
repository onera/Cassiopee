#=============================================================================
# Apps requires:
# ELSAPROD variable defined in environment
# CASSIOPEE
#=============================================================================
import os
from setuptools import setup

prod = os.getenv("ELSAPROD") or "xx"

# setup ======================================================================
setup(
    name=("Cassiopee-" if os.getenv("CASSIOPEE_DIST_PREFIX") else "") + "Apps",
    version="4.2",
    description="Application modules",
    author="ONERA",
    url="https://onera.github.io/Cassiopee/",
    packages=['Apps', 'Apps.Chimera', 'Apps.Fast', 'Apps.Mesh', 'Apps.Coda'],
    package_dir={"":"."}
)
