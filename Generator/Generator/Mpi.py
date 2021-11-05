# parallel Generator
from . import PyTree
import Converter.Mpi4py as Cmpi

def bbox(t):
    """Return the bounding box of a pytree."""
    bb = PyTree.bbox(t)
    bb1 = Cmpi.allreduce(bb, Cmpi.MIN)
    bb2 = Cmpi.allreduce(bb, Cmpi.MAX)
    return [bb1[0],bb1[1],bb1[2],bb2[3],bb2[4],bb2[5]]
