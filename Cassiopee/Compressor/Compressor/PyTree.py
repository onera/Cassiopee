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
    """Write interpolation coefficients for unsteady computations."""
    Compressor.writeUnsteadyCoefs(iteration, indices, filename, loc, format)

# Remplace les coordonnees d'une grille cartesienne par un noeud CartesianData
#
# if layers not None, only communicate the desired number of layers
# bboxDict is dict with the zones of t as keys and their specific bboxes as key values, used when layers not None
# if subr, the tree subregions are kept during the exchange
def compressCartesian(t, bbox=[], layers=None, subr=True, tol=1.e-10):
    """For Cartesian grids, replace Grid Coordinates with a compressed node."""
    tp = Internal.copyRef(t)
    _compressCartesian(tp, bbox=bbox, layers=layers, subr=subr, tol=tol)
    return tp

# compress a gridCoordinates container
def _compressCartesian__(z, ztype, gc, tol=1.e-10):
    xp = Internal.getNodeFromName1(gc, 'CoordinateX')
    yp = Internal.getNodeFromName1(gc, 'CoordinateY')
    zp = Internal.getNodeFromName1(gc, 'CoordinateZ')
    if xp is None: return False
    if yp is None: return False
    if zp is None: return False
    if Internal.getNodeFromName1(xp, 'ZData'): return False # already compress by something else
    if Internal.getNodeFromName1(yp, 'ZData'): return False
    if Internal.getNodeFromName1(zp, 'ZData'): return False
    xp = xp[1].ravel(order='K')
    yp = yp[1].ravel(order='K')
    zp = zp[1].ravel(order='K')
    ni = ztype[1]; nj = ztype[2]; nk = ztype[3]
    x0 = xp[0]; y0 = yp[0]; z0 = zp[0]
    if ni > 3: hi = xp[3]-xp[2]
    elif ni > 2: hi = xp[2]-xp[1]
    elif ni > 1: hi = xp[1]-x0
    else: hi = 1.
    cartesian = True
    if ni > 1:
        if abs(xp[1] - xp[0] - hi) > tol: cartesian = False
        if abs(yp[1] - y0) > tol: cartesian = False
        if abs(zp[1] - z0) > tol: cartesian = False
    if ni > 2:
        if abs(xp[2] - xp[1] - hi) > tol: cartesian = False
        if abs(yp[1] - y0) > tol: cartesian = False
        if abs(zp[1] - z0) > tol: cartesian = False
    if ni > 3: # safe from here
        if abs(xp[3] - xp[2] - hi) > tol: cartesian = False
        if abs(yp[2] - y0) > tol: cartesian = False
        if abs(zp[2] - z0) > tol: cartesian = False
    if ni > 4:
        if abs(xp[4] - xp[3] - hi) > tol: cartesian = False
        if abs(yp[3] - y0) > tol: cartesian = False
        if abs(zp[3] - z0) > tol: cartesian = False
    if ni > 5:
        if abs(xp[5] - xp[4] - hi) > tol: cartesian = False
        if abs(yp[4] - y0) > tol: cartesian = False
        if abs(zp[4] - z0) > tol: cartesian = False

    if nj > 3: hj = yp[3*ni]-yp[2*ni]
    elif nj > 2: hj = yp[2*ni]-yp[ni]
    elif nj > 1: hj = yp[ni]-y0
    else: hj = 1.
    if nj > 1:
        if abs(yp[ni] - yp[0] - hj) > tol: cartesian = False
        if abs(xp[ni] - x0) > tol: cartesian = False
        if abs(zp[ni] - z0) > tol: cartesian = False
    if nj > 2:
        if abs(yp[2*ni] - yp[ni] - hj) > tol: cartesian = False
        if abs(xp[ni] - x0) > tol: cartesian = False
        if abs(zp[ni] - z0) > tol: cartesian = False
    if nj > 3:
        if abs(yp[3*ni] - yp[2*ni] - hj) > tol: cartesian = False
        if abs(xp[2*ni] - x0) > tol: cartesian = False
        if abs(zp[2*ni] - z0) > tol: cartesian = False
    if nj > 4:
        if abs(yp[4*ni] - yp[3*ni] - hj) > tol: cartesian = False
        if abs(xp[3*ni] - x0) > tol: cartesian = False
        if abs(zp[3*ni] - z0) > tol: cartesian = False
    if nj > 5:
        if abs(yp[5*ni] - yp[4*ni] - hj) > tol: cartesian = False
        if abs(xp[4*ni] - x0) > tol: cartesian = False
        if abs(zp[4*ni] - z0) > tol: cartesian = False

    if nk > 3: hk = zp[3*ni*nj]-zp[2*ni*nj]
    elif nk > 2: hk = zp[2*ni*nj]-zp[ni*nj]
    elif nk > 1: hk = zp[ni*nj]-z0
    else: hk = 1.
    if nk > 1:
        if abs(zp[ni*nj] - zp[0] - hk) > tol: cartesian = False
        if abs(xp[ni*nj] - x0) > tol: cartesian = False
        if abs(yp[ni*nj] - y0) > tol: cartesian = False
    if nk > 2:
        if abs(zp[2*ni*nj] - zp[ni*nj] - hk) > tol: cartesian = False
        if abs(xp[ni*nj] - x0) > tol: cartesian = False
        if abs(yp[ni*nj] - y0) > tol: cartesian = False
    if nk > 3:
        if abs(zp[3*ni*nj] - zp[2*ni*nj] - hk) > tol: cartesian = False
        if abs(xp[2*ni*nj] - x0) > tol: cartesian = False
        if abs(yp[2*ni*nj] - y0) > tol: cartesian = False
    if nk > 4:
        if abs(zp[4*ni*nj] - zp[3*ni*nj] - hk) > tol: cartesian = False
        if abs(xp[3*ni*nj] - x0) > tol: cartesian = False
        if abs(yp[3*ni*nj] - y0) > tol: cartesian = False
    if nk > 5:
        if abs(zp[5*ni*nj] - zp[4*ni*nj] - hk) > tol: cartesian = False
        if abs(xp[4*ni*nj] - x0) > tol: cartesian = False
        if abs(yp[4*ni*nj] - y0) > tol: cartesian = False

    #print(cartesian, abs(zp[2*ni*nj] - zp[ni*nj] - hk), abs(zp[3*ni*nj] - zp[2*ni*nj] - hk))
    if cartesian:
        px = Internal.createUniqueChild(gc, 'CoordinateX', 'DataArray_t', value=[0.]*10) # important for skeleton read
        Internal.createChild(px, 'ZData', 'DataArray_t', value=[6., x0, hi, float(ni), float(nj), float(nk)])
        py = Internal.createUniqueChild(gc, 'CoordinateY', 'DataArray_t', value=[0.]*10)
        Internal.createChild(py, 'ZData', 'DataArray_t', value=[6., y0, hj, float(ni), float(nj), float(nk)])
        pz = Internal.createUniqueChild(gc, 'CoordinateZ', 'DataArray_t', value=[0.]*10)
        Internal.createChild(pz, 'ZData', 'DataArray_t', value=[6., z0, hk, float(ni), float(nj), float(nk)])
        cd = Internal.createChild(gc, 'CartesianData', 'DataArray_t', value=[x0,y0,z0,hi,hj,hk])
        Internal.createChild(cd, 'ZData', 'DataArray_t', value=[x0,y0,z0,hi,hj,hk]) # to avoid recompression of CartesianData

    return cartesian

# Si la zone est cartesienne :
# Ajoute un noeud CartesianData a la zone
# remplace Coordinates par des noeuds avec des champs de taille 10
def _compressCartesian(t, bbox=[], layers=None, subr=True, tol=1.e-10):
    """For Cartesian grids, replace Grid Coordinates with a compressed node."""
    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getZoneDim(z)
        if ztype[0] == 'Unstructured': continue
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if gc is None: continue
        cartesian = _compressCartesian__(z, ztype, gc, tol)
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
        if gc is not None: _compressCartesian__(z, ztype, gc, tol)

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
    if c is None: return None # no CartesianData
    c = Internal.getValue(c)
    if c is None: return None # suppose skeleton zone
    x0 = c[0]; y0 = c[1]; z0 = c[2]
    hi = c[3]; hj = c[4]; hk = c[5]
    tmp = G.cart((x0,y0,z0), (hi,hj,hk), (ztype[1], ztype[2], ztype[3]))
    gct = Internal.getNodeFromName1(tmp, Internal.__GridCoordinates__)
    if gc is not None:
        gct[0] = gc[0]
        cn = Internal.getNodeFromName1(gc, 'CoordinateX')
        if cn is not None and cn[1] is None: return None # suppose skeleton zone
    if gc is None: Internal._addChild(z, gct, pos=0)
    else:
        Internal._rmNodesFromName1(z, gc[0])
        Internal._addChild(z, gct, pos=0)
    # cartesianData is automatically removed
    return None

# Si la zone n'est pas skeleton et contient un noeud CartesianData :
# Reconstruit CoordinateX, CoordinateY, CoordinateZ
# Supprime CartesianData
def _uncompressCartesian(t):
    """For Cartesian grids, recreate Grid Coordinates from compressed zones."""
    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getZoneDim(z)
        if ztype[0] == 'Unstructured': continue
        gc = Internal.getNodeFromName1(z, 'GridCoordinates#Init')
        if gc is not None: _uncompressCartesian__(z, ztype, gc)
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if gc is not None: _uncompressCartesian__(z, ztype, gc)
    return None

# Node data compression
# ctype=0: compress with sz
# ctype=1: compress with zfp
# ctype=2: compress cellN (lossless)
# ctype=3: compress basic element connectivity (lossless)
# ctype=4: compress ngon connectivity (losless)
# ctype=5: compress with fpc (lossless)
# ctype=6: reserve pour compressCartesian (lossless)
def _packNode(node, tol=1.e-8, ctype=0):
    if Internal.getNodeFromName1(node, 'ZData') is not None: return None # already compressed node
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
    elif ctype == 5: # fpc
        #print('compress', node[0], node[1].shape, flush=True)
        ret = Compressor.compressor.compressFpc(node[1])
        shape = [5.,0.,float(ret[2])]+list(ret[0])
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
        elif ctype == 5: # fpc
            #print('uncompress', node[0], shape, flush=True)
            ret = Compressor.compressor.uncompressFpc((shape,node[1],iscorder))
            node[1] = ret
            Internal._rmNodesFromName1(node, 'ZData')
        else:
            raise ValueError("unpackNode: unknown compression type.")
    return None

# compressCoords of zones
def _compressCoords(t, tol=1.e-8, ctype=0):
    """Compress coordinates lossless or with a relative tolerance."""
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
    """Compress coordinates lossless or with a relative tolerance."""
    tp = Internal.copyRef(t)
    _compressCoords(tp, tol, ctype)
    return tp

def _compressFields(t, tol=1.e-8, ctype=0, varNames=None):
    """Compress fields lossless or with a relative tolerance."""
    zones = Internal.getZones(t)
    for z in zones:
        if varNames is None:
            # Compress all FlowSolution_t containers
            FS = Internal.getNodesFromType1(z, 'FlowSolution_t')
            fields = []
            for c in FS: fields += Internal.getNodesFromType1(c, 'DataArray_t')
            for f in fields:
                if Internal.getNodeFromName1(f, 'ZData') is None:
                    _packNode(f, tol, ctype)
        else:
            # Compress only given variables
            for v in varNames:
                varname = v.split(':', 1)
                container = Internal.__FlowSolutionNodes__
                if len(varname) == 2 and varname[0] == 'centers':
                    container = Internal.__FlowSolutionCenters__
                    varname = varname[1]
                elif len(varname) == 2 and varname[0] == 'nodes':
                    varname = varname[1]
                else: varname = v
                FS = Internal.getNodeFromName1(z, container)
                f = Internal.getNodeFromName1(FS, varname)
                if f is not None:
                    if Internal.getNodeFromName1(f, 'ZData') is None:
                        _packNode(f, tol, ctype)
                else: print("Warning: compressFields: field %s not found."%v)
    return None

def compressFields(t, tol=1.e-8, ctype=0, varNames=None):
    """Compress fields lossless or with a relative tolerance."""
    tp = Internal.copyRef(t)
    _compressFields(tp, tol, ctype, varNames)
    return tp

# Compresse un cellN 0,1,2
def _compressCellN(t, varNames=['cellN']):
    """Compress cellN (0,1,2) lossless on 2 bits."""
    zones = Internal.getZones(t)
    for name in varNames:
        spl = name.split(':')
        if len(spl) == 2: spl = spl[1]
        else: spl = name
        for z in zones:
            cellNs = Internal.getNodesFromName2(z, spl)
            for cellN in cellNs: _packNode(cellN, 0., 2)
    return None

def compressCellN(t, varNames=['cellN']):
    """Compress cellN (0,1,2) lossless on 2 bits."""
    tp = Internal.copyRef(t)
    _compressCellN(tp, varNames)
    return tp

# Compress Elements_t (elts basiques ou NGONs)
def _compressElements(t):
    """Compress lossless Element connectivities."""
    zones = Internal.getZones(t)
    for z in zones:
        elts = Internal.getNodesFromType1(z, 'Elements_t')
        for e in elts:
            eltno = e[1][0]
            if eltno == 22 or eltno == 23: # NGON
                n = Internal.getNodeFromName1(e, 'ElementConnectivity')
                _packNode(n, 0, 4)
            else:
                _, net = Internal.eltNo2EltName(eltno)
                n = Internal.getNodeFromName1(e, 'ElementConnectivity')
                _packNode(n, net, 3)
    return None

def compressElements(t):
    """Compress lossless Element connectivity."""
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
        for c in GC+FS:
            fields += Internal.getNodesFromType1(c, 'DataArray_t')
        for f in fields:
            _unpackNode(f)
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

# compresse le plus possible en lossless (sauf le cartesien)
def _compressAll(t):
    """Compress coords, fields and connectivity (lossless)."""
    _compressCellN(t)
    _compressCoords(t, ctype=5)
    _compressFields(t, ctype=5)
    _compressElements(t)
    return None

def compressAll(t):
    """Compress coords, fields and connectivity (lossless)."""
    tp = Internal.copyRef(t)
    _compressAll(tp)
    return tp