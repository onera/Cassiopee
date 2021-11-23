# parallel Generator
from . import PyTree
import Converter.Mpi as Cmpi
import numpy

def bbox(t):
    """Return the bounding box of a pytree."""
    bb = PyTree.bbox(t)
    bb = numpy.array(bb, dtype=numpy.float64)
    bb1 = numpy.empty((bb.size), dtype=numpy.float64)
    bb2 = numpy.empty((bb.size), dtype=numpy.float64)
    Cmpi.Allreduce(bb, bb1, Cmpi.MIN)
    Cmpi.Allreduce(bb, bb2, Cmpi.MAX)
    return [bb1[0],bb1[1],bb1[2],bb2[3],bb2[4],bb2[5]]
