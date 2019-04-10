#
# Python Interface to compress solutions from PyTrees
#
from . import Compressor
import numpy
import Converter.Internal as Internal
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
    newPeriodicity=interpData[2]
    coefs=interpData[3]
    oldRcvIndices = ref[0]
    oldDonorIndices = ref[1]
    oldPeriodicity = ref[2]
    oldCoefs=ref[3]
    if loc == 'face':
        newFaceDir = interpData[4]
        oldFaceDir = ref[4]
    else:
        newFaceDir = None
    r1 = numpy.in1d(newRcvIndices,oldRcvIndices)
    r2 = numpy.in1d(newRcvIndices,oldRcvIndices, invert=True)
    r3 = numpy.in1d(oldRcvIndices,newRcvIndices, invert=True)
    alreadyExistingIndices = newRcvIndices[r1] # indices in oldRcvIndices and newRcvIndices
    newIndices = newRcvIndices[r2]             # indices in newRcvIndices but not in oldRcvIndices
    indicesToDelete = oldRcvIndices[r3]        # indices in oldRcvIndices but not in newRcvIndices
    flag = {} # key: rcvIndex, value : Flag
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

def container__(flag,newPos,indDonor,periodicity,coefs,faceDir):
    if flag == 0 and faceDir == None: return [flag,(int)(indDonor),(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]
    elif flag == 0: return [flag,(int)(indDonor),(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 1 and faceDir == None: return [flag]+[(float)(c) for c in coefs[newPos]]
    elif flag == 1: return [flag]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 2 and faceDir == None: return [flag,(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]
    elif flag == 2: return [flag,(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 3 and faceDir == None: return [flag,(int)(indDonor)]+[(float)(c) for c in coefs[newPos]]
    elif flag == 3: return [flag,(int)(indDonor)]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])]
    elif flag == 5: return [flag]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])] # only for 'face'
    elif flag == 6: return [flag,(int)(periodicity[newPos])]+[(float)(c) for c in coefs[newPos]]+[(int)(faceDir[newPos])] # only for 'face'

#------------------------------------------------------------------------------
# writeUnsteadyCoefs
#------------------------------------------------------------------------------
def writeUnsteadyCoefs(iteration, indices, filename, loc,format="b"):
    """write interpolation coefficients for unsteady computation."""
    Compressor.writeUnsteadyCoefs(iteration, indices, filename,loc,format)
    
# Remplace les coordonnees d'une grille cartesienne par un noeud UserDefined
def compressCartesian(t):
    tp = Internal.copyRef(t)
    _compressCartesian(tp)
    return tp

def _compressCartesian(t):
    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getZoneDim(z)
        if ztype[0] == 'Unstructured': continue
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        if gc is None: continue
        xp = Internal.getNodeFromName1(gc, 'CoordinateX')
        yp = Internal.getNodeFromName1(gc, 'CoordinateY')
        zp = Internal.getNodeFromName1(gc, 'CoordinateZ')
        if xp is None: continue
        if yp is None: continue
        if zp is None: continue
        xp = xp[1].ravel(order='K')
        yp = yp[1].ravel(order='K')
        zp = zp[1].ravel(order='K')
        ni = ztype[1]; nj = ztype[2]; nk = ztype[3]
        x0 = xp[0]; y0 = yp[0]; z0 = zp[0]
        hi = xp[1]-x0; 
        if nj > 1: hj = yp[ni]-y0
        else: hj = 1.
        if nk > 1: hk = zp[ni*nj]-z0
        else: hk = 1.
        Internal._rmNodesFromName(gc, 'CoordinateX')
        Internal._rmNodesFromName(gc, 'CoordinateY')
        Internal._rmNodesFromName(gc, 'CoordinateZ')
        Internal.createChild(z, 'CartesianData', 'DataArray_t', value=[x0,y0,z0,hi,hj,hk])
    return None
    
# uncompress Cartesian
def uncompressCartesian(t):
    tp = Internal.copyRef(t)
    _uncompressCartesian(tp)
    return tp
    
def _uncompressCartesian(t):
    import Generator.PyTree as G
    zones = Internal.getZones(t)
    for z in zones:
        ztype = Internal.getZoneDim(z)
        if ztype[0] == 'Unstructured': continue
        c = Internal.getNodeFromName1(z, 'CartesianData')
        if c is None: continue
        c = Internal.getValue(c)
        x0 = c[0]; y0 = c[1]; z0 = c[2]
        hi = c[3]; hj = c[4]; hk = c[5]
        tmp = G.cart((x0,y0,z0), (hi,hj,hk), (ztype[1], ztype[2], ztype[3]))
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        gct = Internal.getNodeFromName1(tmp, Internal.__GridCoordinates__) 
        if gc is None: Internal._addChild(z, gct)
        else: 
            Internal._rmNodesFromName(z, Internal.__GridCoordinates__)
            Internal._addChild(z, gct)
    return None
        
