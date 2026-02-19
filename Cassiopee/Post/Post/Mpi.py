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

    fcenters, fareas = G.getFaceCentersAndAreas(t)
    centers = G.getCellCenters(t, fcenters, fareas)
    zones = Internal.getZones(t)
    for i, zone in enumerate(zones):
        cc = centers[i]
        if cc is None: continue
        fsolc = Internal.getNodeFromName(zone, Internal.__FlowSolutionCenters__)
        Internal.createNode('CCx', 'DataArray_t', cc[0], None, fsolc)
        Internal.createNode('CCy', 'DataArray_t', cc[1], None, fsolc)
        Internal.createNode('CCz', 'DataArray_t', cc[2], None, fsolc)

    # exchange fields and cell centers
    allNames = fldNames
    allNames.append('CCx')
    allNames.append('CCy')
    allNames.append('CCz')
    rfields = X.exchangeFields(t, fldNames)

    # get comm list
    zgc = Internal.getNodeFromType1(zone, 'ZoneGridConnectivity_t')
    if zgc is None: raise ValueError('ZoneGridConnectivity not found.')
    comms = Internal.getNodesFromType1(zgc, 'GridConnectivity1to1_t')
    if comms is None: raise ValueError('GridConnectivity1to1 not found.')
    ptlists = []
    for comm in comms:
        ptlist = Internal.getNodeFromName(comm, 'PointList')[1]
        ptlists.append(ptlist)

    # compute lsq gradients
    parRun = 1
    t = P.computeGradLSQ(t, fldNames, parRun, fcenters, ptlists, rfields)
    # TODO(Imad): delete cell centers from tree
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
# Parallel computeGrad2 for NGON zones
# BCMatch must be set in t
#=============================================================================
def _computeGrad2(t, var, withCellN=True):
    """Compute the gradient of a variable defined in array."""
    if Cmpi.size == 1: return P._computeGrad2(t, var, withCellN)

    import Converter.converter
    import Post

    vare = var.split(':')
    if len(vare) > 1: vare = vare[1]

    # Compute graph of match
    procDict = Cmpi.getProcDict(t)
    graph = Cmpi.computeGraph(t, type='match', procDict=procDict)
    zones = Internal.getZones(t)
    export = {}

    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Unstructured' and dim[3] == 'NGON':
            # adaptation needed by actual computeGrad2
            Internal._adaptNGon42NGon3(z)
            Internal._adaptNFace2PE(z, remove=False)

            # get face values
            GCs = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for gc in GCs:
                donor = Internal.getValue(gc)
                PL = Internal.getBCFaceNode(z, gc)[1] # PointList
                PLD = Internal.getBCFaceNode(z, gc, donor=True)[1] # PointListDonor
                fld = Converter.converter.extractBCMatchNG(z, PL, [vare],
                                                           Internal.__GridCoordinates__,
                                                           Internal.__FlowSolutionNodes__,
                                                           Internal.__FlowSolutionCenters__)
                oppNode = procDict[donor]
                n = [donor, z[0], fld, PLD.ravel('k')]
                if oppNode not in export: export[oppNode] = [n]
                else: export[oppNode] += [n]
        elif dim[0] == 'Structured':
            fields = C.getField(vare, z, api=3)
            # get face values
            GCs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for gc in GCs:
                donor = Internal.getValue(gc)
                #PL = Internal.getBCFaceNode(z, gc)[1] # PointRange
                #PLD = Internal.getBCFaceNode(z, gc, donor=True)[1] # PointRangeDonor
                prr = Internal.getNodeFromName1(gc, 'PointRange')
                prd = Internal.getNodeFromName1(gc, 'PointRangeDonor')
                wr = Internal.range2Window(prr[1])
                wd = Internal.range2Window(prd[1])
                (iminR,imaxR,jminR,jmaxR,kminR,kmaxR) = wr
                (iminD,imaxD,jminD,jmaxD,kminD,kmaxD) = wd
                niR = dim[1]; njR = dim[2]; nkR = dim[3]
                tri = Internal.getNodeFromName1(gc, 'Transform')
                tri = Internal.getValue(tri)
                (t1,t2,t3) = tri
                [indR,fld]  = Converter.converter.extractBCMatchStruct(fields,(iminD,jminD,kminD,imaxD,jmaxD,kmaxD),
                                                                       (iminR,jminR,kminR,imaxR,jmaxR,kmaxR),
                                                                       (niR,njR,nkR),(t1,t2,t3))
                oppNode = procDict[donor]
                n = [donor, z[0], fld, PLD.ravel('k')]
                if oppNode not in export: export[oppNode] = [n]
                else: export[oppNode] += [n]

    # sendrecv
    recvDatas = Cmpi.sendRecv(export, graph)
    #print(Cmpi.rank, 'recv', recvDatas)

    # Mean on faces (we must find the opposite face from donor name)
    indices = {}; BCField = {}
    for i in recvDatas:
        for n in recvDatas[i]:
            # donor is supposed to have a unique matching match
            (donor, source, fld, PLD) = n
            z = Internal.getNodeFromName2(t, donor)
            zn = z[0]
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Unstructured' and dim[3] == 'NGON':
                fld1 = Converter.converter.buildBCMatchFieldNG(z, PLD, fld, [vare],
                                                               Internal.__GridCoordinates__,
                                                               Internal.__FlowSolutionNodes__,
                                                               Internal.__FlowSolutionCenters__)
            elif dim[0] == 'Structured':
                fields = C.getAllFields(z, 'centers', api=3)
                fld1 = Converter.converter.buildBCMatchFieldStruct(fields, PLD, fld, None)
            else:
                raise(TypeError, "computeGrad2: invalid grid type.")
            #print(Cmpi.rank, 'PLD', PLD, indices, flush=True)
            if zn not in indices: indices[zn] = PLD
            else: indices[zn] = numpy.concatenate((indices[zn], PLD))
            if zn not in BCField: BCField[zn] = fld1[1][0].ravel('k')
            else: BCField[zn] = numpy.concatenate((BCField[zn], fld1[1][0].ravel('k')))

            # debug extract
            #if Cmpi.rank == 0:
            #    zp1 = T.subzone(z, PL.ravel('k'), type='faces')
            #    C.convertPyTree2File(zp1, 'out1.cgns')
            #G._getVolumeMap(z)

    for z in zones:
        zn = z[0]
        # Test if vol present
        vol = None
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if cont is not None:
            vol = Internal.getNodeFromName1(cont, 'vol')
            if vol is not None: vol = vol[1]

        # Test if cellN present
        cellN = None
        if withCellN:
            if cont is not None:
                cellN  = Internal.getNodeFromName1(cont, 'cellN')
                if cellN is not None: cellN = cellN[1]
        else: cellN = None

        # get BCDataSet
        isghost = Internal.getNodeFromType1(z, 'Rind_t')
        if isghost is None: # not a ghost cells zone : add BCDataSet
            zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
            if zoneBC is not None:
                BCs = Internal.getNodesFromType1(zoneBC, 'BC_t')
                for b in BCs:
                    datas = Internal.getBCDataSet(z, b)
                    inds = Internal.getBCFaceNode(z, b)
                    if datas != [] and inds != []:
                        bcf = None
                        for i in datas:
                            if i[0] == vare: bcf = i; break
                        if bcf is not None:
                            indsp = inds[1].ravel(order='K')
                            bcfp = bcf[1].ravel(order='K')
                            if zn not in indices: indices[zn] = indsp
                            else: indices[zn] = numpy.concatenate((indices[zn], indsp))
                            if zn not in BCField: BCField[zn] = bcfp
                            else: BCField[zn] = numpy.concatenate((BCField[zn], bcfp))

        f = C.getField(var, z, api=1)[0]
        x = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]

        if f != []:
            if zn in indices: inds = indices[zn]
            else: inds = None
            if zn in BCField: bcf = BCField[zn]
            else: bcf = None

            centers = Post.computeGrad2(x, f, vol, cellN, indices=inds, BCField=bcf)
            C.setFields([centers], z, 'centers')

    return None

#=============================================================================
# Parallel computeFiv2 for NGON zones
# BCMatch must be set in t
#=============================================================================
def _computeDiv2(t, var, withCellN=True):
    """Compute the divergence of a variable defined in array."""
    if Cmpi.size == 1: return P._computeDiv2(t, var, withCellN)

    import Converter.converter
    import Post

    vare = var.split(':')
    if len(vare) > 1: vare = vare[1]

    # Compute graph of match
    procDict = Cmpi.getProcDict(t)
    graph = Cmpi.computeGraph(t, type='match', procDict=procDict)
    zones = Internal.getZones(t)
    export = {}

    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Unstructured' and dim[3] == 'NGON':
            # adaptation needed by actual computeDiv2
            Internal._adaptNGon42NGon3(z)
            Internal._adaptNFace2PE(z, remove=False)

            # get face values
            GCs = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for gc in GCs:
                donor = Internal.getValue(gc)
                PL = Internal.getBCFaceNode(z, gc)[1] # PointList
                PLD = Internal.getBCFaceNode(z, gc, donor=True)[1] # PointListDonor
                fld = Converter.converter.extractBCMatchNG(z, PL, [vare],
                                                           Internal.__GridCoordinates__,
                                                           Internal.__FlowSolutionNodes__,
                                                           Internal.__FlowSolutionCenters__)
                oppNode = procDict[donor]
                n = [donor, z[0], fld, PLD.ravel('k')]
                if oppNode not in export: export[oppNode] = [n]
                else: export[oppNode] += [n]
        elif dim[0] == 'Structured':
            fields = C.getField(vare, z, api=3)
            # get face values
            GCs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            for gc in GCs:
                donor = Internal.getValue(gc)
                #PL = Internal.getBCFaceNode(z, gc)[1] # PointRange
                #PLD = Internal.getBCFaceNode(z, gc, donor=True)[1] # PointRangeDonor
                prr = Internal.getNodeFromName1(gc, 'PointRange')
                prd = Internal.getNodeFromName1(gc, 'PointRangeDonor')
                wr = Internal.range2Window(prr[1])
                wd = Internal.range2Window(prd[1])
                (iminR,imaxR,jminR,jmaxR,kminR,kmaxR) = wr
                (iminD,imaxD,jminD,jmaxD,kminD,kmaxD) = wd
                niR = dim[1]; njR = dim[2]; nkR = dim[3]
                tri = Internal.getNodeFromName1(gc, 'Transform')
                tri = Internal.getValue(tri)
                (t1,t2,t3) = tri
                [indR,fld] = Converter.converter.extractBCMatchStruct(fields,(iminD,jminD,kminD,imaxD,jmaxD,kmaxD),
                                                                      (iminR,jminR,kminR,imaxR,jmaxR,kmaxR),
                                                                      (niR,njR,nkR),(t1,t2,t3))
                oppNode = procDict[donor]
                n = [donor, z[0], fld, PLD.ravel('k')]
                if oppNode not in export: export[oppNode] = [n]
                else: export[oppNode] += [n]

    # sendrecv
    recvDatas = Cmpi.sendRecv(export, graph)

    # Mean on faces (we must find the opposite face from donor name)
    indices = {}; BCField = {}
    for i in recvDatas:
        for n in recvDatas[i]:
            # donor is supposed to have a unique matching match
            (donor, source, fld, PLD) = n
            z = Internal.getNodeFromName2(t, donor)
            zn = z[0]
            dim = Internal.getZoneDim(z)
            if dim[0] == 'Unstructured' and dim[3] == 'NGON':
                fld1 = Converter.converter.buildBCMatchFieldNG(z, PLD, fld, [vare],
                                                               Internal.__GridCoordinates__,
                                                               Internal.__FlowSolutionNodes__,
                                                               Internal.__FlowSolutionCenters__)
            elif dim[0] == 'Structured':
                fields = C.getAllFields(z, 'centers', api=3)
                fld1 = Converter.converter.buildBCMatchFieldStruct(fields, PLD, fld, None)
            else:
                raise(TypeError, "computeDiv2: invalid grid type.")
            if zn not in indices: indices[zn] = PLD
            else: indices[zn] = numpy.concatenate((indices[zn], PLD))
            if zn not in BCField: BCField[zn] = fld1[1][0].ravel('k')
            else: BCField[zn] = numpy.concatenate((BCField[zn], fld1[1][0].ravel('k')))

    for z in zones:
        zn = z[0]
        # Test if vol present
        vol = None
        cont = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
        if cont is not None:
            vol = Internal.getNodeFromName1(cont, 'vol')
            if vol is not None: vol = vol[1]

        # Test if cellN present
        cellN = None
        if withCellN:
            if cont is not None:
                cellN  = Internal.getNodeFromName1(cont, 'cellN')
                if cellN is not None: cellN = cellN[1]
        else: cellN = None

        # get BCDataSet
        isghost = Internal.getNodeFromType1(z, 'Rind_t')
        if isghost is None: # not a ghost cells zone : add BCDataSet
            zoneBC = Internal.getNodesFromType1(z, 'ZoneBC_t')
            if zoneBC is not None:
                BCs = Internal.getNodesFromType1(zoneBC, 'BC_t')
                for b in BCs:
                    datas = Internal.getBCDataSet(z, b)
                    inds = Internal.getBCFaceNode(z, b)
                    if datas != [] and inds != []:
                        bcf = None
                        for i in datas:
                            if i[0] == vare: bcf = i; break
                        if bcf is not None:
                            indsp = inds[1].ravel(order='K')
                            bcfp = bcf[1].ravel(order='K')
                            if zn not in indices: indices[zn] = indsp
                            else: indices[zn] = numpy.concatenate((indices[zn], indsp))
                            if zn not in BCField: BCField[zn] = bcfp
                            else: BCField[zn] = numpy.concatenate((BCField[zn], bcfp))

        f = C.getField(var, z, api=1)[0]
        x = C.getFields(Internal.__GridCoordinates__, z, api=1)[0]

        if f != []:
            if zn in indices: inds = indices[zn]
            else: inds = None
            if zn in BCField: bcf = BCField[zn]
            else: bcf = None

            centers = Post.computeDiv2(x, f, vol, cellN, indices=inds, BCField=bcf)
            C.setFields([centers], z, 'centers')

    return None
