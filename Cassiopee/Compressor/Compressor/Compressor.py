"""Module for CFD solution compression.
"""
__version__ = '4.1'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud, Xavier Juvigny"

from . import compressor
import numpy
import Converter.Internal as Internal

#==============================================================================
# deltaIndex
#==============================================================================
def deltaIndex(index, ref):
    """Return the delta between index and ref."""
    try: r1 = numpy.isin(index, ref)
    except: r1 = numpy.is1d(index, ref)
    r1 = r1.astype(Internal.E_NpyInt)
    try: r2 = numpy.isin(ref, index)
    except: r2 = numpy.is1d(ref, index)
    r2 = r2.astype(Internal.E_NpyInt)
    return compressor.deltaIndex(index, ref, r1, r2)

#==============================================================================
# writeUnsteadyCoefs
#==============================================================================
def writeUnsteadyCoefs(iteration, indices, filename, loc, format="b"):
    """Write interpolation coefficients for unsteady computation."""
    compressor.writeUnsteadyCoefs(iteration, indices, filename, loc, format)
