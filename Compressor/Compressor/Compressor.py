"""Module for CFD solution compression.
"""
__version__ = '2.8'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud"

import compressor
import numpy

#==============================================================================
# deltaIndex
#==============================================================================
def deltaIndex(index, ref):
    """Return the delta between index and ref."""
    r1 = numpy.in1d(index, ref)
    r1 = r1.astype(numpy.int32)
    r2 = numpy.in1d(ref, index)
    r2 = r2.astype(numpy.int32)
    return compressor.deltaIndex(index, ref, r1, r2)

#==============================================================================
# writeUnsteadyCoefs
#==============================================================================
def writeUnsteadyCoefs(iteration, indices, filename, loc, format="b"):
    """Write interpolation coefficients for unsteady computation."""
    compressor.writeUnsteadyCoefs(iteration, indices, filename, loc, format)

#==============================================================================
# Serialize/compress
# method=0: pickle
#==============================================================================
def pack(a, method=0):
    if method == 0:
        import cPickle as pickle# best for now
        return pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        import cPickle as pickle # best for now
        return pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL)
   
#==============================================================================
# Unserialize/decompress
#==============================================================================
def unpack(a, method=0):
    if method == 0:
        import cPickle as pickle# best for now
        return pickle.loads(a)
    else:
        import cPickle as pickle # best for now
        return pickle.loads(a)
