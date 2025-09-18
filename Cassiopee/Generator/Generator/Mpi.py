# parallel Generator
from . import PyTree
from . import Generator
import numpy
import Converter.Mpi as Cmpi
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G
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

def getRegularityMap(t, addGC=False):
    """Return the regularity map in an array.
    Usage: getRegularityMap(t)"""
    if addGC: t = Cmpi.addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    t = C.TZGC1(t, 'centers', True, Generator.getRegularityMap)
    if addGC: t = Internal.rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return t

def _getRegularityMap(t, addGC=False):
    """Return the regularity map in an array."""
    if addGC: Cmpi._addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    C._TZGC1(t, 'centers', False, Generator.getRegularityMap)
    if addGC: Internal._rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return None

def getAngleRegularityMap(t, addGC=False):
    """Return the regularity map in an array (wrt angles).
    Usage: getAngleRegularityMap(t)"""
    if addGC: t = Cmpi.addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    t = C.TZGC1(t, 'centers', True, Generator.getAngleRegularityMap)
    if addGC: t = Internal.rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return t

def _getAngleRegularityMap(t, addGC=False):
    """Return the regularity map in an array (wrt angles)."""
    if addGC: Cmpi._addGhostCells(t, t, 1, adaptBCs=0, modified=[], fillCorner=1)
    C._TZGC1(t, 'centers', False, Generator.getAngleRegularityMap)
    if addGC: Internal._rmGhostCells(t, t, 1, adaptBCs=0, modified=[])
    return None

#========================================================
# check mesh
#========================================================
def getMeshFieldInfo__(m, field, critValue, verbose):
    fmin  = 1.e32; fsum  = 0; fmax  = -1.; fcrit = 0; size  = 0
    info = 'INFO %s: min = %1.2e, max = %1.2e, mean = %1.2e, crit(%s %s %s) = %s cells out of %s | %2.2f%% (%s)'

    dictInfo = {}

    for z in Internal.getZones(m):
        f = Internal.getNodeFromName(z, field)[1]

        zname = z[0]

        size_loc  = numpy.size(f)
        fcrit_loc = numpy.count_nonzero(f<critValue) if field == 'vol' else numpy.count_nonzero(f>critValue)
        fmin_loc  = numpy.min(f)
        fmax_loc  = numpy.max(f)
        fsum_loc  = numpy.sum(f)

        fmin   = min(fmin_loc, fmin)
        fmax   = max(fmax_loc, fmax)
        fsum  += fsum_loc
        fcrit += fcrit_loc
        size  += size_loc

        dictInfo[zname] = [Cmpi.rank, fmin_loc, fmax_loc, fsum_loc, fcrit_loc, size_loc]

    Cmpi.barrier()

    dictInfo = Cmpi.gather(dictInfo, root=0)

    if Cmpi.rank == 0:
        dictInfo = {k:v for d in dictInfo for k,v in d.items()}
        for zname in dictInfo:
            rank_loc, fmin_loc, fmax_loc, fsum_loc, fcrit_loc, size_loc = dictInfo[zname]
            if verbose == 2 or (verbose == 1 and fcrit_loc > 0):
                print(info%(field.upper(),fmin_loc,fmax_loc,fsum_loc/float(size_loc),field,'<' if field == 'vol' else '>',critValue,fcrit_loc,size_loc,fcrit_loc/float(size_loc)*100,"rank %s - %s"%(rank_loc,zname)))

    fmin  = Cmpi.allreduce(fmin,  op=Cmpi.MIN)
    fmax  = Cmpi.allreduce(fmax,  op=Cmpi.MAX)
    fsum  = Cmpi.allreduce(fsum,  op=Cmpi.SUM)
    fcrit = Cmpi.allreduce(fcrit, op=Cmpi.SUM)
    size  = Cmpi.allreduce(size,  op=Cmpi.SUM)

    fmean = fsum/float(size)

    if Cmpi.rank == 0 and (verbose == 2 or (verbose == 1 and fcrit_loc > 0)):
        print('#'*(len(field)+7))
        print(info%(field.upper(),fmin,fmax,fsum/float(size),field,'<' if field == 'vol' else '>',critValue,fcrit,size,fcrit/float(size)*100,'GLOBAL'))
        print('#'*(len(field)+7)+'\n')

    Cmpi.barrier()

    return fmin, fmax, fmean, fcrit

def checkMesh(m, critVol=0., critOrtho=15., critReg=0.1, critAngReg=15., addGC=False, verbose=0):
    """Return information on mesh quality."""

    Cmpi.barrier()

    G._getVolumeMap(m)
    vmin,vmax,vmean,vcrit = getMeshFieldInfo__(m, 'vol', critVol, verbose)
    Internal._rmNodesFromName(m, 'vol')

    G._getOrthogonalityMap(m)
    omin,omax,omean,ocrit = getMeshFieldInfo__(m, 'orthogonality', critOrtho, verbose)
    Internal._rmNodesFromName(m, 'orthogonality')

    _getRegularityMap(m, addGC)
    rmin,rmax,rmean,rcrit = getMeshFieldInfo__(m, 'regularity', critReg, verbose)
    Internal._rmNodesFromName(m, 'regularity')

    _getAngleRegularityMap(m, addGC)
    amin,amax,amean,acrit = getMeshFieldInfo__(m, 'regularityAngle', critAngReg, verbose)
    Internal._rmNodesFromName(m, 'regularityAngle')

    return {'vmin':vmin,'vmax':vmax,'vmean':vmean,'vcrit':vcrit,
            'rmin':rmin,'rmax':rmax,'rmean':rmean,'rcrit':rcrit,
            'amin':amin,'amax':amax,'amean':amean,'acrit':acrit,
            'omin':omin,'omax':omax,'omean':omean,'ocrit':ocrit}

#=============================================================================================
# mesh adaptation
#=============================================================================================
def adaptMesh(a, indicator="indicator", hook=None, dim=3, conformize=False,
              splitInfos=None):

    if splitInfos is None:
        print("Warning: no input provided to adapt the parts of the mesh. They will be adapted independently.", flush=True)
    a = Internal.getZones(a)[0]
    dimZ = Internal.getZoneDim(a)
    eltType=dimZ[3]
    if eltType != 'NGON':
        print("adaptMesh: input mesh must be NGON v4. No adaptation performed.", flush=True)
        return a

    return G.adaptMesh__(a, indicator=indicator, hook=hook, dim=dim,
                         conformize=conformize, splitInfos=splitInfos)
