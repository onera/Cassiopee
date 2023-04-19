# parallel Generator
from . import PyTree
from . import Generator
import numpy
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
import numpy

from mpi4py import MPI

def bbox(t):
    """Return the bounding box of a pytree."""
    bb = PyTree.bbox(t)
    bb = numpy.array(bb, dtype=numpy.float64)
    bb1 = numpy.empty((bb.size), dtype=numpy.float64)
    bb2 = numpy.empty((bb.size), dtype=numpy.float64)
    Cmpi.Allreduce(bb, bb1, Cmpi.MIN)
    Cmpi.Allreduce(bb, bb2, Cmpi.MAX)
    return [bb1[0],bb1[1],bb1[2],bb2[3],bb2[4],bb2[5]]

def getRegularityMap(t, addGC=False):
    """Return the regularity map in an array.
    Usage: getRegularityMap(t)"""
    if addGC: t = Cmpi.addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    t = C.TZGC(t, 'centers', Generator.getRegularityMap)
    if addGC: t = Internal.rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return t

def _getRegularityMap(t, addGC=False):
    if addGC: Cmpi._addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    C._TZGC(t, 'centers', Generator.getRegularityMap)
    if addGC: Internal._rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return None
    
def getAngleRegularityMap(t, addGC=False):
    """Return the regularity map in an array (wrt angles).
    Usage: getAngleRegularityMap(t)"""
    if addGC: t = Cmpi.addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    t = C.TZGC(t, 'centers', Generator.getAngleRegularityMap)
    if addGC: t = Internal.rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return t

def _getAngleRegularityMap(t, addGC=False):
    if addGC: Cmpi._addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    C._TZGC(t, 'centers', Generator.getAngleRegularityMap)
    if addGC: Internal._rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return None

#========================================================
# Mesh quality informations
#========================================================
def getMeshFieldInfo(m, field, critDict, verbose=True):
    fmin  = 1.e32
    fsum  = 0
    fmax  = -1.
    fcrit = 0
    size  = 0

    for z in Internal.getZones(m):
        f = Internal.getNodeFromName(z, field)[1]

        size_loc  = numpy.size(f)
        fcrit_loc = numpy.count_nonzero(f<critDict[field]) if field == 'vol' else numpy.count_nonzero(f>critDict[field]) 
        fmin_loc  = numpy.min(f)
        fmax_loc  = numpy.max(f)
        fsum_loc  = numpy.sum(f)

        fmin   = min(fmin_loc, fmin)
        fmax   = max(fmax_loc, fmax)
        fsum  += fsum_loc
        fcrit += fcrit_loc
        size  += size_loc

        if verbose:
            print('INFO {} min = {:1.2e}, max = {:1.2e}, mean = {:1.2e}, crit({} {} {}) = {} cells out of {} | {:2.2f}% (rank {} - {})'.format(field.upper()+' : ',fmin_loc,fmax_loc,fsum_loc/float(size_loc),field,'<' if field == 'vol' else '>',critDict[field],fcrit_loc,size_loc,fcrit_loc/float(size_loc)*100,Cmpi.rank,z[0]), flush=True)

    Cmpi.barrier()

    fmin  = Cmpi.allreduce(fmin,  op=Cmpi.MIN)
    fmax  = Cmpi.allreduce(fmax,  op=Cmpi.MAX)
    fsum  = Cmpi.allreduce(fsum,  op=Cmpi.SUM)
    fcrit = Cmpi.allreduce(fcrit, op=Cmpi.SUM)
    size  = Cmpi.allreduce(size,  op=Cmpi.SUM)

    fmean = fsum/float(size)

    if Cmpi.rank == 0 and verbose:
        print('#'*(len(field)+7))
        print('INFO {} min = {:1.2e}, max = {:1.2e}, mean = {:1.2e}, crit({} {} {}) = {} cells out of {} | {:2.2f}% ({})'.format(field.upper()+' : ',fmin,fmax,fmean,field,'<' if field == 'vol' else '>',critDict[field],fcrit,size,fcrit/float(size)*100,'GLOBAL'), flush=True)
        print('#'*(len(field)+7)+'\n')

    Cmpi.barrier()

    return fmin, fmax, fmean, fcrit

def checkMesh(m, critDict={'vol':0., 'regularity':0.1, 'orthogonality':15., 'regularityAngle':15.}, addGC=False, verbose=True):
    """Return information on mesh quality."""

    Cmpi.barrier()

    G._getVolumeMap(m)
    vmin,vmax,vmean,vcrit = getMeshFieldInfo(m, 'vol', critDict, verbose)

    _getRegularityMap(m, addGC)
    rmin,rmax,rmean,rcrit = getMeshFieldInfo(m, 'regularity', critDict, verbose)

    _getAngleRegularityMap(m, addGC)
    amin,amax,amean,acrit = getMeshFieldInfo(m, 'regularityAngle', critDict, verbose)
    
    G._getOrthogonalityMap(m)
    omin,omax,omean,ocrit = getMeshFieldInfo(m, 'orthogonality', critDict, verbose)

    Cmpi.convertPyTree2File(m, 'checkmesh.cgns')

    return {'vmin':vmin,'vmax':vmax,'vmean':vmean,'vcrit':vcrit,
            'rmin':rmin,'rmax':rmax,'rmean':rmean,'rcrit':rcrit,
            'amin':amin,'amax':amax,'amean':amean,'acrit':acrit,
            'omin':omin,'omax':omax,'omean':omean,'ocrit':ocrit}
