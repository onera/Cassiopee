"""Module for CFD solution compression.
"""
__version__ = '2.9'
__author__ = "Stephanie Peron, Christophe Benoit, Pascal Raud"

from . import compressor
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
    """Serialize or compress a."""
    if method == 0:
        try: import cPickle as pickle # best for now
        except: import pickle 
        return pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL)
    elif method == 1:
        try: import cPickle as pickle # best for now
        except: import picle
        import zlib
        return zlib.compress(pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL), level=1)
    else:
        try: import cPickle as pickle # best for now
        except: import pickle
        return pickle.dumps(a, protocol=pickle.HIGHEST_PROTOCOL)
   
#==============================================================================
# Unserialize/decompress
#==============================================================================
def unpack(a, method=0):
    """Deserialize or decompress a."""
    if method == 0:
        try: import cPickle as pickle # best for now
        except: import pickle
        return pickle.loads(a)
    elif method == 1:
        try: import cPickle as pickle # best for now
        except: import pickle
        import zlib
        return pickle.loads(zlib.decompress(a))
    else:
        try: import cPickle as pickle # best for now
        except: import pickle
        return pickle.loads(a)
