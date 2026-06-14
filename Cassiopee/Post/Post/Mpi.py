# Interface pour MPI

import Converter.Mpi as Cmpi
from . import PyTree as P
import Converter.Internal as Internal
import Converter.PyTree as C
import KCore.Vector as Vector
import numpy
import XCore.PyTree as X
import Generator.PyTree as G

#==============================================================================
# Compute grad with LSQ method
#==============================================================================
def computeGradLSQ(t, fldNames):
    if Cmpi.size == 1: return P.computeGradLSQ(t, fldNames)

    if isinstance(fldNames, str):
        fldNames = [fldNames]
    fldNames = [fldName.split(':')[-1] for fldName in fldNames]

    fcenters, fareas = G.getFaceCentersAndAreas(t)
    centers = G.getCellCenters(t, fcenters, fareas)
    zones = Internal.getZones(t)
    for i, zone in enumerate(zones):
        cc = centers[i]
        if cc is None: continue
        fsolc = Internal.getNodeFromName(zone, Internal.__FlowSolutionCenters__)
        if fsolc is None: raise ValueError("computeGradLSQ: FlowSolutionCenters not found.")
        Internal.newDataArray('CCx', cc[0], fsolc)
        Internal.newDataArray('CCy', cc[1], fsolc)
        Internal.newDataArray('CCz', cc[2], fsolc)

    # exchange fields and cell centers
    allFldNames = fldNames + ['CCx', 'CCy', 'CCz']
    rfields = X.exchangeFields(t, allFldNames)

    # get comm list
    zgc = Internal.getNodeFromType1(zone, 'ZoneGridConnectivity_t')
    if zgc is None: raise ValueError('computeGradLSQ: ZoneGridConnectivity not found.')
    comms = Internal.getNodesFromType1(zgc, 'GridConnectivity1to1_t')
    if comms is None: raise ValueError('computeGradLSQ: GridConnectivity1to1 not found.')
    ptlists = []
    for comm in comms:
        ptlist = Internal.getNodeFromName(comm, 'PointList')[1]
        ptlists.append(ptlist)

    # compute lsq gradients
    parRun = 1
    t = P.computeGradLSQ(t, fldNames, parRun, fcenters, ptlists, rfields)
    return t

#==============================================================================
# extractMesh
# IN: t: maillage source distribue
# IN: extractionMesh: maillage de destination distribue
# IN: graph: graph d'intersection si deja calcule
#==============================================================================
def extractMesh(t, extractionMesh, order=2, extrapOrder=1,
                constraint=40., tol=1.e-6, mode='robust', hook=None, graph=None):
    if Cmpi.size == 1: return P.extractMesh(t, extractionMesh, order, extrapOrder,
                                            constraint, tol, mode, hook)
    if graph is None:
        tb = Cmpi.createBBoxTree(t)
        tb2 = Cmpi.createBBoxTree(extractionMesh)
        graph = Cmpi.computeGraph(tb, type='bbox3', t2=tb2)
    tl = Cmpi.addXZones(t, graph)
    tl = Cmpi.convert2PartialTree(tl)
    ext = Cmpi.convert2PartialTree(extractionMesh)
    # print info
    nztl = len(Internal.getZones(tl))
    nzext = len(Internal.getZones(ext))
    print('Rank %d has %d source zones and %d destination zones.'%(Cmpi.rank, nztl, nzext))
    ext = P.extractMesh(tl, ext, order=order, extrapOrder=extrapOrder, constraint=constraint, tol=tol, mode=mode,
                        hook=hook)
    return ext

#==============================================================================
def integ(t, var=''):
    """Integral of fields defined in t."""
    if t is not None: ret = P.integ(t, var)
    else: ret = 0.
    ret  = numpy.array(ret, dtype=numpy.float64)
    ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    Cmpi.Allreduce(ret, ret1, Cmpi.SUM)
    return ret1.tolist()

#==============================================================================
def integNorm(t, var=''):
    """Integral of fields times normal."""
    if t is not None: ret = P.integNorm(t, var)
    else: ret = 0.
    ret  = numpy.array(ret, dtype=numpy.float64)
    ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    Cmpi.Allreduce(ret, ret1, Cmpi.SUM)
    return [ret1.tolist()]

#==============================================================================
def integNormProduct(t, vector=[]):
    """Integral of fields product normal."""
    if t is not None: ret = P.integNormProduct(t, vector)
    else: ret = 0.
    ret = numpy.array(ret, dtype=numpy.float64)
    ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    Cmpi.Allreduce(ret, ret1, Cmpi.SUM)
    return ret1.tolist()

#==============================================================================
def integMoment(t, center=(0.,0.,0.), vector=[]):
    """Integral of moments."""
    if t is not None: ret = P.integMoment(t, center, vector)
    else: ret = 0.
    ret = numpy.array(ret, dtype=numpy.float64)
    ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    Cmpi.Allreduce(ret, ret1, Cmpi.SUM)
    return ret1.tolist()

#==============================================================================
def integMomentNorm(t, center=(0.,0.,0.), var=''):
    """Integral of moments."""
    if t is not None: ret = P.integMomentNorm(t, center, var)
    else: ret = 0.
    ret = numpy.array(ret, dtype=numpy.float64)
    ret1 = numpy.empty(ret.shape, dtype=numpy.float64)
    Cmpi.Allreduce(ret, ret1, Cmpi.SUM)
    return ret1.tolist()

#=============================================================================
# Parallel streamline2 : dans la direction de l'ecoulement uniquement (dir=1)
#=============================================================================
def streamLine2(t, X0, vector, N=2000, eps=1.e-2, maxCompt=20):
    """Compute a streamline starting from (x0,y0,z0) given
    a list of arrays containing 'vector' information."""

    if Cmpi.size == 1: return P.streamLine2(t, X0, vector, N, eps, maxCompt)
    out = []; compt = 0
    while len(X0) > 0 and compt < maxCompt:
        ret = P.streamLine2(t, X0, vector, N=N, dir=1, eps=eps)
        for z in ret: z[0] = z[0]+'_%d'%Cmpi.rank

        # Get new pool (supprime les streamlines degenerees)
        X0 = []; ret2 = []
        for z in ret:
            P0 = C.getValue(z, 'GridCoordinates', -1)
            P1 = C.getValue(z, 'GridCoordinates', 1)
            dP = Vector.sub(P0, P1)
            l = Vector.norm2(dP)
            if l >= 1.e-10:
                Pts = P0 # last point
                X0.append(tuple(Pts))
                ret2.append(z)
        #print('>> New pool', X0)
        out += ret2

        # Communicate and merge pool
        b = Cmpi.allgather(X0)
        X0 = []
        for i in b: X0 += i
        #print('>> New pool after com', X0)
        print('it=%d pool length=%d'%(compt,len(X0)))
        compt += 1

    return out

#=============================================================================
# Helper functions for computeGrad2 and computeDiv2
#=============================================================================
def _getVolAndCellN__(z, vol, cellN, withCellN=True):
    fsolc = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    if fsolc is not None:
        vol = Internal.getNodeFromName1(fsolc, 'vol')
        if vol is not None: vol = vol[1]
        cellN = Internal.getNodeFromName1(fsolc, 'cellN')
        if cellN is not None: cellN = cellN[1]
    return None

def _addBCDataSet2ExchangeBCMatchData__(z, varList, indices, BCField):
    # Add BCDataSet data if zone is not a ghost cell zone
    isghost = Internal.getNodeFromType1(z, 'Rind_t')
    zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
    if (isghost is not None) or zoneBC is None: return

    zn = z[0]
    nvars = len(varList)
    BCs = Internal.getNodesFromType1(zoneBC, 'BC_t')
    for b in BCs:
        datas = Internal.getBCDataSet(z, b)
        inds = Internal.getBCFaceNode(z, b)
        if datas == [] or inds == []: continue
        bcf = [None for _ in range(nvars)]
        for i in datas:
            for v in range(nvars):
                if i[0] == varList[v] and bcf[v] is None: bcf[v] = i
        if all(bcfv is not None for bcfv in bcf):
            indsp = inds[1].ravel(order='K')
            bcfp = [bcfv[1].ravel(order='K') for bcfv in bcf]
            if zn not in indices: indices[zn] = indsp
            else: indices[zn] = numpy.concatenate((indices[zn], indsp))
            for v, bcfpv in enumerate(bcfp):
                BCFieldv = BCField[v]
                if zn not in BCFieldv: BCFieldv[zn] = bcfpv
                else: BCFieldv[zn] = numpy.concatenate((BCFieldv[zn], bcfpv))
    return None

#=============================================================================
# Parallel computeGrad2 for STRUCT and NGON zones
# BCMatch must be set in t
#=============================================================================
def _computeGrad2(t, var, withCellN=True):
    """Compute the gradient of a variable defined in array."""
    if Cmpi.size == 1:
        return P._computeGrad2(t, var, withCellN=withCellN)

    import Connector.Mpi as Xmpi
    import Post

    # Exchange BCMatch data
    varList = [var.split(':')[-1]]
    indices, BCField = Xmpi.exchangeBCMatchData(t, varList)

    zones = Internal.getZones(t)
    for z in zones:
        zn = z[0]
        f = C.getField(var, z, api=3)[0]
        if f == []: continue

        # Test if vol and cellN are present
        vol, cellN = None, None
        _getVolAndCellN__(z, vol, cellN, withCellN=withCellN)

        # Add BCDataSet data to exchange BCMatch data
        _addBCDataSet2ExchangeBCMatchData__(z, varList, indices, BCField)

        x = C.getFields(Internal.__GridCoordinates__, z, api=3)[0]
        inds = indices[zn] if zn in indices else None
        bcf = BCField[0][zn] if zn in BCField[0] else None

        gradFlds = Post.computeGrad2(x, f, vol, cellN, indices=inds, BCField=bcf)
        C.setFields([gradFlds], z, 'centers')

    return None

#=============================================================================
# Parallel computeDiv2 for STRUCT and NGON zones
# BCMatch must be set in t
#=============================================================================
def _computeDiv2(t, var, rmVar=False):
    """Compute the divergence of a variable defined in array."""
    if Cmpi.size == 1:
        return P._computeDiv2(t, var, rmVar=rmVar)

    import Connector.Mpi as Xmpi
    import Post

    # Exchange BCMatch data
    varList = [var.split(':')[-1] + d for d in ["X", "Y", "Z"]]
    indices, BCField = Xmpi.exchangeBCMatchData(t, varList)

    zones = Internal.getZones(t)
    for z in zones:
        zn = z[0]
        f = C.getFields("centers", z, vars=varList, api=3)[0]
        if f == []: continue

        # Test if vol and cellN are present
        vol, cellN = None, None
        _getVolAndCellN__(z, vol, cellN)

        # Add BCDataSet data to exchange BCMatch data
        _addBCDataSet2ExchangeBCMatchData__(z, varList, indices, BCField)

        x = C.getFields(Internal.__GridCoordinates__, z, api=3)[0]
        BCFieldX = BCField[0]; BCFieldY = BCField[1]; BCFieldZ = BCField[2]
        inds = indices[zn] if zn in indices else None
        bcfx = BCFieldX[zn] if zn in BCFieldX else None
        bcfy = BCFieldY[zn] if zn in BCFieldY else None
        bcfz = BCFieldZ[zn] if zn in BCFieldZ else None

        divFld = Post.computeDiv2(
            x, f, vol, cellN, indices=inds,
            BCFieldX=bcfx, BCFieldY=bcfy, BCFieldZ=bcfz
        )
        C.setFields([divFld], z, 'centers')

    # Conditional clean up of partial derivatives
    if rmVar:
        varWLocList = [var + d for d in ["X", "Y", "Z"]]
        C._rmVars(t, varWLocList)
        C._rmBCDataVars(t, varList)

    return None
