#=============================================================================
# Apps requires:
# ELSAPROD variable defined in environment
# CASSIOPEE
#=============================================================================
import os
from setuptools import setup
import KCore.Dist as Dist

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

# Write build info to file
Dist.writeBuildInfo("Apps/App.py")

# setup ======================================================================
setup()
