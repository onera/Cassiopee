"""Module for CFD solution compression.
"""
from . import Compressor
import numpy
import Converter.Internal as Internal
import Converter.PyTree as C
__version__ = Compressor.__version__

#------------------------------------------------------------------------------
# deltaInterpolations
#------------------------------------------------------------------------------
def deltaInterpolations(interpData, ref, loc='cell'):
    """Return the delta between index and ref."""
    # Build Flag
    #  0: new point (full storage)
    # -1: deleted point (storage of rcv index)
    #  1: existing point with other donor (storage of donor index and coefficients)
    #  2: existing point with other coefficients (storage of coefficients)
    newRcvIndices = interpData[0]
    newDonorIndices = interpData[1]
    newPeriodicity = interpData[2]
    coefs = interpData[3]
    oldRcvIndices = ref[0]
    oldDonorIndices = ref[1]
    oldPeriodicity = ref[2]
    oldCoefs = ref[3]
    if loc == 'face':
        newFaceDir = interpData[4]
        oldFaceDir = ref[4]
    else:
        newFaceDir = None
    r1 = numpy.in1d(newRcvIndices, oldRcvIndices)
    r2 = numpy.in1d(newRcvIndices, oldRcvIndices, invert=True)
    r3 = numpy.in1d(oldRcvIndices, newRcvIndices, invert=True)
    alreadyExistingIndices = newRcvIndices[r1] # indices in oldRcvIndices and newRcvIndices
    newIndices = newRcvIndices[r2]             # indices in newRcvIndices but not in oldRcvIndices
    indicesToDelete = oldRcvIndices[r3]        # indices in oldRcvIndices but not in newRcvIndices
    #flag = {} # key: rcvIndex, value : Flag
    storedInterpData = {} # key : rcvIndex, value : storage flag + interpolation data to write
    # Deals with new indices:
    for rcvIndex in newIndices:
        newPos = numpy.where(newRcvIndices==rcvIndex)[0][0]
        indDonor = newDonorIndices[newPos]
        storageFlag = 0 # new point
        storedInterpData[(int)(rcvIndex)] = container__(0,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
    # Deals with indices to delete
    for rcvIndex in indicesToDelete:
        storageFlag = -1 # deleted point
        storedInterpData[(int)(rcvIndex)]=[storageFlag]
    # Deals with already existing indices
    for rcvIndex in alreadyExistingIndices:
        newPos = numpy.where(newRcvIndices==rcvIndex)[0][0]
        indDonor = newDonorIndices[newPos]
        oldPos = numpy.where(oldRcvIndices==rcvIndex)[0][0]
        oldIndDonor = oldDonorIndices[oldPos]
        if indDonor == oldIndDonor: # same donor index
            if newPeriodicity[newPos] == oldPeriodicity[oldPos]:
                if loc == 'cell':
                    if numpy.array_equal(coefs[newPos],oldCoefs[oldPos]) == False:
                        storageFlag = 1
                        storedInterpData[(int)(rcvIndex)]= container__(1,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
                    else: # DBG
                        storageFlag = 4 # DBG
                        storedInterpData[(int)(rcvIndex)]=[storageFlag] # DBG
                else:
                    if newFaceDir[newPos] == oldFaceDir[oldPos]:
                        if numpy.array_equal(coefs[newPos],oldCoefs[oldPos]) == False:
                            storageFlag = 1
                            storedInterpData[(int)(rcvIndex)]= container__(1,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
                        else: # DBG
                            storageFlag = 4 # DBG
                            storedInterpData[(int)(rcvIndex)]=[storageFlag] # DBG
                    else:
                        storageFlag = 5
                        storedInterpData[(int)(rcvIndex)]= container__(5,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
            else:
                if loc == 'cell':
                    storageFlag = 2
                    storedInterpData[(int)(rcvIndex)]= container__(2,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
                else:
                    if newFaceDir[newPos] == oldFaceDir[oldPos]:
                        storageFlag = 2
                        storedInterpData[(int)(rcvIndex)]= container__(2,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
                    else:
                        storageFlag = 6
                        storedInterpData[(int)(rcvIndex)]= container__(6,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
        else: # different donor
            if newPeriodicity[newPos] == oldPeriodicity[oldPos]:
                storageFlag = 3
                storedInterpData[(int)(rcvIndex)] = container__(3,newPos,indDonor,newPeriodicity,coefs,newFaceDir)
            else:
                storageFlag = 0 # full storage as a new point
                storedInterpData[(int)(rcvIndex)] = container__(0,newPos,indDonor,newPeriodicity,coefs,newFaceDir)

    return storedInterpData

def container__(flag, newPos, indDonor, periodicity, coefs, faceDir):
    if flag == 0 and faceDir is None: return [flag,(int)(indDonor),(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]
    elif flag == 0: return [flag,(int)(indDonor),(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 1 and faceDir is None: return [flag]+[(float)(c) for c in coefs[newPos]]
    elif flag == 1: return [flag]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 2 and faceDir is None: return [flag,(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]
    elif flag == 2: return [flag,(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 3 and faceDir is None: return [flag,(int)(indDonor)]+[(float)(c) for c in coefs[newPos]]
    elif flag == 3: return [flag,(int)(indDonor)]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 5: return [flag]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])] # only for 'face'
    elif flag == 6: return [flag,(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])] # only for 'face'

#------------------------------------------------------------------------------
# writeUnsteadyCoefs
#------------------------------------------------------------------------------
def writeUnsteadyCoefs(iteration, indices, filename, loc, format="b"):
    """write interpolation coefficients for unsteady computations."""
    Compressor.writeUnsteadyCoefs(iteration, indices, filename, loc, format)
    
# Remplace les coordonnees d'une grille cartesienne par un noeud CartesianData
#
# if layers not None, only communicate the desired number of layers
# bboxDict is dict with the zones of t as keys and their specific bboxes as key values, used when layers not None
# if subr, the tree subregions are kept during the exchange 
def compressCartesian(t, bbox=[], layers=None, subr=True):
    """For Cartesian grids, replace Grid Coordinates with a UserDefined CartesianData node."""
    tp = Internal.copyRef(t)
    _compressCartesian(tp, bbox=bbox, layers=layers, subr=subr)
    return tp

# compress a gridCorddinates container
def _compressCartesian__(z, ztype, gc):
    xp = Internal.getNodeFromName1(gc, 'CoordinateX')
    yp = Internal.getNodeFromName1(gc, 'CoordinateY')
    zp = Internal.getNodeFromName1(gc, 'CoordinateZ')
    if xp is None: return False
    if yp is None: return False
    if zp is None: return False
    xp = xp[1].ravel(order='K')
    yp = yp[1].ravel(order='K')
    zp = zp[1].ravel(order='K')
    ni = ztype[1]; nj = ztype[2]; nk = ztype[3]
    x0 = xp[0]; y0 = yp[0]; z0 = zp[0]
    hi = xp[1]-x0
    cartesian = True
    if ni > 2:
        if abs(xp[2] - xp[1] - hi) > 1.e-10: cartesian = False
        if abs(yp[1] - y0) > 1.e-10: cartesian = False
        if abs(zp[1] - z0) > 1.e-10: cartesian = False
    if ni > 3:
        if abs(xp[3] - xp[2] - hi) > 1.e-10: cartesian = False
        if abs(yp[3] - y0) > 1.e-10: cartesian = False
        if abs(zp[3] - z0) > 1.e-10: cartesian = False

    if nj > 1: hj = yp[ni]-y0
    else: hj = 1.
    if nj > 2:
        if abs(yp[2*ni] - yp[ni] - hj) > 1.e-10: cartesian = False
        if abs(xp[ni] - x0) > 1.e-10: cartesian = False
        if abs(zp[ni] - z0) > 1.e-10: cartesian = False
    if nj > 3:
        if abs(yp[3*ni] - yp[2*ni] - hj) > 1.e-10: cartesian = False
        if abs(xp[2*ni] - x0) > 1.e-10: cartesian = False
        if abs(zp[2*ni] - z0) > 1.e-10: cartesian = False
        
    if nk > 1: hk = zp[ni*nj]-z0
    else: hk = 1.
    if nk > 2: 
        if abs(zp[2*ni*nj] - zp[ni*nj] - hk) > 1.e-10: cartesian = False
        if abs(xp[ni*nj] - x0) > 1.e-10: cartesian = False
        if abs(yp[ni*nj] - y0) > 1.e-10: cartesian = False
    if nk > 3: 
        if abs(zp[3*ni*nj] - zp[2*ni*nj] - hk) > 1.e-10: cartesian = False
        if abs(xp[2*ni*nj] - x0) > 1.e-10: cartesian = False
        if abs(yp[2*ni*nj] - y0) > 1.e-10: cartesian = False

    if cartesian:
        Internal.createUniqueChild(gc, 'CoordinateX', 'DataArray_t', value=[0.]*10) # important for skeleton read
        Internal.createUniqueChild(gc, 'CoordinateY', 'DataArray_t', value=[0.]*10)
        Internal.createUniqueChild(gc, 'CoordinateZ', 'DataArray_t', value=[0.]*10)
        Internal.createChild(gc, 'CartesianData', 'DataArray_t', value=[x0,y0,z0,hi,hj,hk])
        
    return cartesian

# Si la zone est cartesienne :
# Ajoute un noeud CartesianData a la zone
# remplace Coordinates par des noeuds avec des champs de taille 10
def _compressCartesian(t, bbox=[], layers=None, subr=True):
    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getZoneDim(z)
        if ztype[0] == 'Unstructured': continue
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if gc is None: continue
        cartesian = _compressCartesian__(z, ztype, gc)
        #print('cartesian?=', cartesian)

        # traitement layers
        if cartesian and layers is not None:
            ni = ztype[1]; nj = ztype[2]; nk = ztype[3]
            xp = Internal.getNodeFromName1(gc, 'CoordinateX')
            yp = Internal.getNodeFromName1(gc, 'CoordinateY')
            zp = Internal.getNodeFromName1(gc, 'CoordinateZ')
            x0 = xp[0]; y0 = yp[0]; z0 = zp[0]
            hi = xp[1]-x0
            if nj > 1: hj = yp[ni]-y0
            else: hj = 1.
            if nk > 1: hk = zp[ni*nj]-z0
            else: hk = 1.

            vars  = C.getVarNames(z, excludeXYZ=True)[0]
            # align bbox with the gridzone
            bbox[0] = numpy.round((bbox[0]-x0)/hi)*hi+x0
            bbox[3] = numpy.round((bbox[3]-x0)/hi)*hi+x0
            # Add the number of layers to the bbox coordinates in every directions and
            # reduce the data volume to an acceptable region (i.e. information cannot be outside the zone)
            x0 = max(xp[0], bbox[0]-layers*hi)
            x1 = min(xp[ni-1], bbox[3]+layers*hi)
            # translate the first and last coordinates of the block to send in terms of a number of layers
            # useful for the numpy slices
            xmin = int(numpy.round((x0 - xp[0])/hi))
            xmax = int(numpy.round((x1- xp[0])/hi))

            bbox[1] = numpy.round((bbox[1]-y0)/hj)*hj+y0
            bbox[4] = numpy.round((bbox[4]-y0)/hj)*hj+y0
            y0 = max(yp[0], bbox[1]-layers*hj)
            y1 = min(yp[ni*nj-1], bbox[4]+layers*hj)
            ymin = int(numpy.round((y0 - yp[0])/hj))
            ymax = int(numpy.round((y1- yp[0])/hj))

            if nk > 1:
                # perform the same treatment in the z-direction
                bbox[2] = numpy.round((bbox[2]-z0)/hk)*hk+z0
                bbox[5] = numpy.round((bbox[5]-z0)/hk)*hk+z0
                z0 = max(zp[0], bbox[2]-layers*hk)
                z1 = min(zp[ni*nj*nk-1], bbox[5]+layers*hk)
                zmin = int(numpy.round((z0 - zp[0])/hk))
                zmax = int(numpy.round((z1- zp[0])/hk))
                for var in vars:
                    if 'centers:' in var:
                        var = var.replace('centers:', '')
                        var = Internal.getNodeFromName(z, var)
                        # order='F' is VERY important
                        Internal.setValue(var, numpy.array(var[1][xmin:xmax,ymin:ymax,zmin:zmax],order='F'))
                    else:
                        # if the data are stored in nodes, we need to go one step further for the numpy slices
                        var = Internal.getNodeFromName(z, var)
                        Internal.setValue(var, numpy.array(var[1][xmin:xmax+1,ymin:ymax+1,zmin:zmax+1],order='F'))
            else:
                for var in vars:
                    if 'centers:' in var:
                        var = var.replace('centers:', '')
                        var = Internal.getNodeFromName(z, var)
                        Internal.setValue(var, numpy.array(var[1][xmin:xmax,ymin:ymax],order='F'))
                    else:
                        var = Internal.getNodeFromName(z, var)
                        Internal.setValue(var, numpy.array(var[1][xmin:xmax+1,ymin:ymax+1],order='F'))

            # Replace the old size information by the new ones
            z[1][0,0] = xmax-xmin+1
            z[1][1,0] = ymax-ymin+1
            if nk > 1: z[1][2,0] = zmax-zmin+1

            if not subr: Internal._rmNodesByType(z, 'ZoneSubRegion_t') # delete subregions (e.g. ID, IBCD in tc)

        # Compress GridInit is present
        gc = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
        if gc is not None: _compressCartesian__(z, ztype, gc)

    return None
    
# uncompress Cartesian
def uncompressCartesian(t):
    """For Cartesian grids, recreate Grid Coordinates from compressed zones."""
    tp = Internal.copyRef(t)
    _uncompressCartesian(tp)
    return tp

def _uncompressCartesian__(z, ztype, gc):
    import Generator.PyTree as G
    c = Internal.getNodeFromName1(gc, 'CartesianData')
    if c is None: return None
    c = Internal.getValue(c)
    x0 = c[0]; y0 = c[1]; z0 = c[2]
    hi = c[3]; hj = c[4]; hk = c[5]
    tmp = G.cart((x0,y0,z0), (hi,hj,hk), (ztype[1], ztype[2], ztype[3]))
    gct = Internal.getNodeFromName1(tmp, Internal.__GridCoordinates__)
    if gc is not None:
        gct[0] = gc[0]
        cn = Internal.getNodeFromName1(gc, 'CoordinateX')
        if cn is not None and cn[1] is None: return None # suppose skeleton zone
    if gc is None: Internal._addChild(z, gct)
    else:
        Internal._rmNodesFromName1(z, gc[0])
        Internal._addChild(z, gct)
    # cartesianData is automatically removed
    return None

# Si la zone n'est pas skeleton et contient un noeud CartesianData :
# Reconstruit CoordinateX, CoordinateY, CoordinateZ
# Supprime CartesianData
def _uncompressCartesian(t):
    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getZoneDim(z)
        if ztype[0] == 'Unstructured': continue
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if gc is not None: _uncompressCartesian__(z, ztype, gc)
        gc = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
        if gc is not None: _uncompressCartesian__(z, ztype, gc)
    
    return None

# ctype=0: compress with sz
# ctype=1: compress with zfp
# ctype=2: compress cellN
def _packNode(node, tol=1.e-8, ctype=0):
    if ctype == 0: # sz
        from . import sz
        ret = sz.pack(node[1], {'relBoundRatio':tol})
        #ret = sz.pack(n, {'absErrBound':tol, 'errorBoundMode':2})
        shape = [0.,tol,float(ret[2])]+list(ret[0])
        node[1] = ret[1]
        Internal._createUniqueChild(node, 'ZData', 'DataArray_t', value=shape)
    elif ctype == 1: # zfp
        from . import zfp
        ret = zfp.pack(node[1], accuracy=tol)
        shape = [1.,tol,float(ret[2])]+list(ret[0])
        node[1] = ret[1]
        Internal._createUniqueChild(node, 'ZData', 'DataArray_t', value=shape)
    elif ctype == 2: # cellN
        ret = Compressor.compressor.compressCellN(node[1])
        iscorder = not numpy.isfortran(node[1])
        shape = [2.,tol,float(iscorder)]+list(ret[0])
        node[1] = ret[1]
        Internal._createUniqueChild(node, 'ZData', 'DataArray_t', value=shape)
    elif ctype == 3: # Elements basiques
        net = int(tol)
        ret = Compressor.compressor.compressIndices((net,node[1]))
        iscorder = not numpy.isfortran(node[1])
        shape = [3.,0.,float(iscorder),net,node[1].size]
        node[1] = ret[2]
        Internal._createUniqueChild(node, 'ZData', 'DataArray_t', value=shape)
    elif ctype == 4: # Elements NGONs
        ret = Compressor.compressor.compressNGonIndices(node[1])
        iscorder = not numpy.isfortran(node[1])
        shape = [4.,0.,float(iscorder),ret[0]]
        node[1] = ret[1]
        Internal._createUniqueChild(node, 'ZData', 'DataArray_t', value=shape)
    else:
        raise ValueError("packNode: unknow compression type.")
    return None

def _unpackNode(node):
    shape = Internal.getNodeFromName1(node, 'ZData')
    if shape is not None and node[1] is not None:
        shape = shape[1]
        ctype = int(shape[0])
        tol = shape[1]
        iscorder = bool(shape[2])
        shape = shape[3:]
        shape = [int(i) for i in shape]
        shape = tuple(shape)
        if ctype == 0: # sz
            from . import sz
            ret = sz.unpack((shape,node[1],iscorder), {})
            node[1] = ret
            Internal._rmNodesFromName1(node, 'ZData')
        elif ctype == 1: # zfp
            from . import zfp
            ret = zfp.unpack((shape,node[1],iscorder), accuracy=tol)
            node[1] = ret
            Internal._rmNodesFromName1(node, 'ZData')
        elif ctype == 2: # cellN
            node[1] = Compressor.compressor.uncompressCellN((shape,node[1],iscorder))
            Internal._rmNodesFromName1(node, 'ZData')
        elif ctype == 3: # Elements (basiques)
            ret = Compressor.compressor.uncompressIndices((shape[0], shape[1], node[1],iscorder))
            node[1] = ret[0]
            Internal._rmNodesFromName1(node, 'ZData')
        elif ctype == 4: # Elements (NGONs)
            ret = Compressor.compressor.uncompressNGonIndices((shape[0],node[1],iscorder))
            node[1] = ret[0]
            Internal._rmNodesFromName1(node, 'ZData')
        else:
            raise ValueError("unpackNode: unknown compression type.")
    return None

# compressFields of zones
def _compressCoords(t, tol=1.e-8, ctype=0):
    """Compress coordinates with a relative tolerance."""
    from . import sz
    zones = Internal.getZones(t)
    for z in zones:
        GC = Internal.getNodesFromType1(z, 'GridCoordinates_t')
        fields = []
        for c in GC: fields += Internal.getNodesFromType1(c, 'DataArray_t')
        for f in fields:
            if Internal.getNodeFromName1(f, 'ZData') is None:
                _packNode(f, tol, ctype)
    return None

def compressCoords(t, tol=1.e-8, ctype=0):
    """Compress coordinates with a relative tolerance."""
    tp = Internal.copyRef(t)
    _compressCoords(tp, tol, ctype)
    return tp

def _compressFields(t, tol=1.e-8, ctype=0):
    """Compress fields with a relative tolerance."""
    zones = Internal.getZones(t)
    for z in zones:
        FS = Internal.getNodesFromType1(z, 'FlowSolution_t')
        fields = []
        for c in FS: fields += Internal.getNodesFromType1(c, 'DataArray_t')
        for f in fields:
            if Internal.getNodeFromName1(f, 'ZData') is None:
                _packNode(f, tol, ctype)
    return None

def compressFields(t, tol=1.e-8, ctype=0):
    """Compress fields with a relative tolerance."""
    tp = Internal.copyRef(t)
    _compressFields(tp, tol, ctype)
    return tp

# Compresse un cellN 0,1,2
def _compressCellN(t):
    """Compress cellN on 2 bits."""
    zones = Internal.getZones(t)
    for z in zones:
        cellNs = Internal.getNodesFromName2(z, 'cellN')
        for cellN in cellNs: _packNode(cellN, 0., 2)
    return None

def compressCellN(t):
    """Compress cellN on 2 bits.""" 
    tp = Internal.copyRef(t)
    _compressCellN(tp)
    return tp

# Compress Elements_t (elts basiques ou NGONs)
def _compressElements(t):
    """Compress Element connectivities."""
    zones = Internal.getZones(t)
    for z in zones:
        elts = Internal.getNodesFromType1(z, 'Elements_t')
        for e in elts:
            eltno = e[1][0]
            if eltno == 22 or eltno == 23: # NGON
                n = Internal.getNodeFromName1(e, 'ElementConnectivity')
                _packNode(n, 0, 4)
            else:                  
                (stype, net) = Internal.eltNo2EltName(eltno)
                n = Internal.getNodeFromName1(e, 'ElementConnectivity')
                _packNode(n, net, 3)
    return None
    
def compressElements(t):
    """Compress Element connectivity.""" 
    tp = Internal.copyRef(t)
    _compressElements(tp)
    return tp

# uncompressFields of zones (si ZData est trouve dans le noeud DataArray_t)
def _uncompressAll(t):
    """Uncompress all compressed data."""
    zones = Internal.getZones(t)
    for z in zones:
        # unpack field
        GC = Internal.getNodesFromType1(z, 'GridCoordinates_t')
        FS = Internal.getNodesFromType1(z, 'FlowSolution_t')
        fields = []
        for c in GC+FS: fields += Internal.getNodesFromType1(c, 'DataArray_t')
        for f in fields: _unpackNode(f)
        # unpack connectivity
        elts = Internal.getNodesFromType1(z, 'Elements_t')
        for e in elts:
            cn = Internal.getNodeFromName1(e, 'ElementConnectivity')
            _unpackNode(cn)
    return None

def uncompressAll(t):
    """Uncompress all compressed data."""
    tp = Internal.copyRef(t)
    _uncompressAll(tp)
    return tp
