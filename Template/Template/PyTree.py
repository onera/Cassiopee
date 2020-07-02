"""Do something cool with pyTrees.
"""
from . import Template
from . import template
import numpy
import Converter.Internal as Internal
__version__ = Template.__version__

# z must be a zone
def pyTreeExample(z):
    return template.pyTreeExample(z, Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__)

# t must be a pyTree
def pyTreeExample1(t):
    return template.pyTreeExample1(t)
