import Converter.PyTree as C
import Transform.PyTree as T
import Converter.Internal as I
import Generator.PyTree as G
import Post.PyTree as P
import Connector.PyTree as X
from . import PyTree as XOR
import time
import numpy
import sys

#==============================================================================
# buildPeriodicCanalOperand

# GOAL : Append the configuration with both neighbor sectors to ensure periodicity of the result 

# IN: t             : 3D NGON mesh

# OUT: returns the tolerance and the max bounding box size
#==============================================================================
def buildPeriodicCanalOperand(zv, angle, JTOL, rotationCenter=[0.,0.,0.], rotationAxis = [1., 0., 0.]):

    #zv = X.connectMatchPeriodic(zv, rotationCenter, angle*rotationAxis, tol=JTOL, dim=3, unitAngle='Degree')

    rotationAxis[0] *= angle ; rotationAxis[1] *= angle ; rotationAxis[2] *= angle
    zvbis = T.rotate(zv, rotationCenter, rotationAxis)

    rotationAxis[0] *= -1. ; rotationAxis[1] *= -1. ; rotationAxis[2] *=  -1.
    zvter = T.rotate(zv, rotationCenter, rotationAxis)

    zv[0] = 'z21'
    zvbis[0] = 'z22'
    zvter[0] = 'z23'
    t2 = C.newPyTree(['Base', zv, zvbis, zvter])
    #t2 = X.connectMatch(t2)
    #C.convertPyTree2File(t2, 'op2.cgns')
    return t2

def transform(z, translation = [0.,0.,0.], rotationCenter=[0.,0.,0.], rotationAxis = [0., 0., 0.]):
    is_trans=False
    if translation[0] != 0. or translation[1] != 0. or translation[2] != 0.:
        is_trans = True;
    is_rot = False
    if rotationAxis[0] != 0. or rotationAxis[1] != 0. or rotationAxis[2] != 0.:
        is_rot = True

    if is_trans == True and is_rot == True:
        print('Input tranform error : translation and rotation are defined both')
    if is_trans == False and is_rot == False:
        print('Input tranform error : no translation nor rotation is defined')

    if is_trans == True: 
        zt = T.translate(z, translation)
    if is_rot == True:
        zt = T.rotate(z, rotationCenter, rotationAxis)
    return zt

# creates a tree with 3 zones : the main and one on each side of it
def build_periodic_operand_2(zi, zone_name, translation=[0.,0.,0.], rotationCenter=[0.,0.,0.], rotationAxis = [0., 0., 0.]):

    z = I.copyRef(zi) # to avoid to modify input zone with connectMatch call

    # 'left zone'
    z2 = transform(z, translation, rotationCenter, rotationAxis)

    z = I.getZones(z)[0]
    z2 = I.getZones(z2)[0]

    z[0]  = zone_name + str(1)
    z2[0] = zone_name + str(2)

    t = C.newPyTree(['Base', z, z2])

    mel_v = XOR.edgeLengthExtrema(z)
    jtol = 0.01*mel_v
    t = X.connectMatch(t, tol=jtol)
    return t

# creates a tree with 3 zones : the main and one on each side of it
def build_periodic_operand_3(zi, zone_name, translation=[0.,0.,0.], rotationCenter=[0.,0.,0.], rotationAxis = [0., 0., 0.]):

    z = I.copyRef(zi) # to avoid to modify input zone with connectMatch call
    # 'left zone'
    z2 = transform(z, translation, rotationCenter, rotationAxis)

    # 'right zone'
    translation[0]  *= -1 ; translation[1]   *= -1 ; translation[2]   *= -1
    rotationAxis[0] *= -1. ; rotationAxis[1] *= -1. ; rotationAxis[2] *= -1.

    z3 = transform(z, translation, rotationCenter, rotationAxis)

    z = I.getZones(z)[0]
    z2 = I.getZones(z2)[0]
    z3 = I.getZones(z3)[0]

    z[0]  = zone_name + str(1)
    z2[0] = zone_name + str(2)
    z3[0] = zone_name + str(3)

    t = C.newPyTree(['Base', z, z2, z3])

    mel_v = XOR.edgeLengthExtrema(z)
    jtol = 0.01*mel_v
    t = X.connectMatch(t, tol=jtol)
    return t

#==============================================================================
# buildPeriodicFeatureOperand

# GOAL : Append the configuration with one neighbor sector to ensure periodicity of the result 

# IN: t             : 3D NGON mesh

# OUT: returns the tolerance and the max bounding box size
#==============================================================================
def buildPeriodicFeatureOperand(zf, angle, JTOL, rotationCenter=[0.,0.,0.], rotationAxis = [1., 0., 0.]):

    rotationAxis[0] *= angle ; rotationAxis[1] *= angle ; rotationAxis[2] *= angle
    zf = X.connectMatchPeriodic(zf, rotationCenter, rotationAxis, tol=JTOL, dim=3, unitAngle='Degree') # to freeze periodic boundaries when agglomerating

    zfbis = T.rotate(zf, rotationCenter, rotationAxis)
    zf[0] = 'z11'
    zfbis[0] = 'z12'
    t1 = C.newPyTree(['Base', zf, zfbis])
    #C.convertPyTree2File(t1, 'op1.cgns')
    return t1

#==============================================================================
# prepareFeature (aperiodic version)

# GOAL : remove singularities on features (overlapping polygons) by agglomeration on the vicinity of the canal's skin to ease the assembly. 

# IN: feature             : 3D NGON mesh of the features to assemble (axisym features)
# IN: canal               : 3D NGON mesh of the canal
# IN: max_overlap_angle   : max angle in radian between two polygons normals to consider them as overlapping
# IN: max_simplify_angle  : sharp edge angle in radian. Under that value, Elligible poylgons will be agglomerated to simplify the cells

# OUT: returns the adapted feature
#==============================================================================
def prepareFeature(feature, canal, max_overlap_angle, max_simplify_angle, treat_externals = 1):
    vf = P.exteriorFaces(canal)
    vf = XOR.convertNGON2DToNGON3D(vf)
    #C.convertPyTree2File(vf, 'vf.cgns')
    #import time
    #t0 = time.time()
    #print('get overlapping faces...')
    res = XOR.getOverlappingFaces(feature, vf, RTOL = 0.15, amax = max_overlap_angle)# rad == 6 deg
    # get pgids for t1 zones only : first par of each pairs
    nb_zones = len(res)
    t1zones_pgids = []
    for i in range(nb_zones):
        t1zones_pgids.append(res[i][0])
    #print('agglomerateCellsWithSpecifiedFaces')
    ag_feature = XOR.agglomerateCellsWithSpecifiedFaces(feature, t1zones_pgids, treat_externals=treat_externals, amax = max_simplify_angle)
    #C.convertPyTree2File(ag_feature, 'ag_feature.cgns')
    return ag_feature

#==============================================================================
# preparePeriodicFeature

# GOAL : remove singularities on features (overlapping polygons) by agglomeration on the vicinity of the canal's skin to ease the assembly. 

# ASSUMPTION : feature and canal have one main zone each. The remaining zones are duplicated and rotated (by rotation_angle) of the main one (for periodicity)
# so we force the agglomeration to be the same for all zones by merging what we get on each zone with getOverlappingFaces 
# and apply the result to all zones when calling agglomerateCellsWithSpecifiedFaces
# this constraint ensure the boolean operation to give a periodic result.

# IN: feature             : 3D NGON mesh of the features to assemble (axisym features)
# IN: canal               : 3D NGON mesh of the canal
# IN: JTOL                : Join tolerance
# IN: rotation_angle      : rotation angle in degree (360 / nb_sectors)
# IN: max_overlap_angle   : max angle in radian between two polygons normals to consider them as overlapping
# IN: max_simplify_angle  : sharp edge angle in radian. Under that value, Elligible poylgons will be agglomerated to simplify the cells
# IN: treat_externals     : Simplification parameter. Does the agglomeration process tend to agglomerate also external faces or not ?

# OUT: returns the adapted feature
#==============================================================================
def preparePeriodicFeature(feature, canal, JTOL, rotation_angle, max_overlap_angle, max_simplify_angle, treat_externals = 1):

    tmp = T.join(canal, tol=JTOL)
    s2 = P.exteriorFaces(tmp)
    s2 = XOR.convertNGON2DToNGON3D(s2)
    #C.convertPyTree2File(s2, 's2.cgns')
    #print('get overlapping faces...')
    res = XOR.getOverlappingFaces(feature, s2, RTOL = 0.15, amax = max_overlap_angle)

    nb_zones = len(res)

    # concatenate the id lists for feature : res[i][0] with i in [0,nb_zones-1]
    idlist = []
    for i in range(nb_zones):
        idlist = numpy.concatenate((idlist, res[i][0]))
    # get rid of duplicate ids
    idlist = list(set(idlist))
    #format
    ids = numpy.empty(len(idlist), I.E_NpyInt)
    ids[:] = idlist[:]  

    # give this list to all zones
    featurezones_pgids = []
    for i in range(nb_zones):
        featurezones_pgids.append(ids)

    #print('agglomerateCellsWithSpecifiedFaces')
    feature = XOR.agglomerateCellsWithSpecifiedFaces(feature, featurezones_pgids, treat_externals=treat_externals, amax = max_simplify_angle) 

    # force periodicity by cloning the most agglomerated zone and replace the less one
    z1s = I.getZones(feature)
    nf1 = XOR.nb_faces(z1s[0])
    nf2 = XOR.nb_faces(z1s[1])

    if nf1 < nf2 : # choose z11
        zb = T.rotate(z1s[0], (0.,0.,0.), (rotation_angle, 0., 0.))
        m = C.getFields(I.__GridCoordinates__, zb)[0]
        C.setFields([m], z1s[1], 'nodes') # replace the mesh in the zone
    else:          # choose z12
        zb = T.rotate(z1s[1], (0.,0.,0.), (-rotation_angle, 0., 0.))
        m = C.getFields(I.__GridCoordinates__, zb)[0]
        C.setFields([m], z1s[0], 'nodes') # replace the mesh in the zone

    return feature


#==============================================================================
# regularizeFeature

# GOAL : remove singularities on features (overlapping polygons) by agglomeration on the vicinity of the canal's skin to ease the assembly. 

# ASSUMPTION : feature and canal have one main zone each. The remaining zones are duplicated and rotated (by rotation_angle) of the main one (for periodicity)
# so we force the agglomeration to be the same for all zones by merging what we get on each zone with getOverlappingFaces 
# and apply the result to all zones when calling agglomerateCellsWithSpecifiedFaces
# this constraint ensure the boolean operation to give a periodic result.

# IN: feature             : 3D NGON mesh of the features to assemble (axisym features)
# IN: canal               : 3D NGON mesh of the canal

# IN: max_overlap_angle   : max angle in radian between two polygons normals to consider them as overlapping
# IN: max_simplify_angle  : sharp edge angle in radian. Under that value, Elligible poylgons will be agglomerated to simplify the cells

# OUT: returns the adapted feature
#==============================================================================
def regularizeFeature(feature, skin, max_overlap_angle, max_simplify_angle):
    res1 = XOR.getOverlappingFaces(feature, skin, RTOL = 0.15, amax = max_overlap_angle)
    res2 = XOR.getCollidingTopFaces(feature, skin, RTOL = 0.15)

    #f = XOR.getFaces(feature, [res2[0]])
    #C.convertPyTree2File(f, 'tops.cgns')

    #f = XOR.getFaces(feature, [res1[0][0]])
    #C.convertPyTree2File(f, 'ovlp.cgns')

    nbz = len(res1) ## == len(res2)
    #print(nbz)
    ids_per_z = []
    for i in range(nbz):
        idlist = []
        if len(res1) > i and len(res1[i]) > 0 : 
            idlist = numpy.concatenate((idlist, res1[i][0]))
        if len(res2) > i : 
            idlist = numpy.concatenate((idlist, res2[i]))
        idlist = list(set(idlist))# get rid of duplicate ids
        ids = numpy.empty(len(idlist), I.E_NpyInt) #format
        ids[:] = idlist[:]
        ids_per_z.append(ids)

    #f = XOR.getFaces(feature, ids_per_z)
    #C.convertPyTree2File(f, 'allf.cgns')

    azt = XOR.agglomerateCellsWithSpecifiedFaces(feature, ids_per_z, treat_externals=1, amax = max_simplify_angle)
    return azt

#==============================================================================
# adaptFirstToSecond(

# GOAL : Makes some checks for compatibility with the solver

# IN: t             : 3D NGON mesh
#
#==============================================================================
def adaptFirstToSecond(comp1, comp2, XVAL, NVAL, Nneigh):
    # XVAL     : levels of subdivision for colliding cells
    # NVAL     : levels of subdivision for the surrounding zone
    # Nneigh   : the surrounding zone is the Nneigh-th neighborhood

    n = XOR.nb_cells(comp1)
    cell_vals = numpy.empty((n,), dtype=I.E_NpyInt)
    cell_vals[:] = 0

    # flag what collides comp2' skin
    s1 = P.exteriorFaces(comp2)
    res = XOR.getCollidingCells(comp1, s1)
    t2xlist = res[0][0]
    #print ('nb of colliding cells : ' + str(len(t2xlist)) + ' over ' + str(n))
    for i in t2xlist:
        cell_vals[i]=XVAL

    # set a value in the neighborhood
    t2neighlist = XOR.getNthNeighborhood(comp1, Nneigh, [t2xlist])
    for i in t2neighlist:
        cell_vals[i]=NVAL#max(NVAL, cell_vals[i])

    comp1 = XOR.adaptCells(comp1, cell_vals, sensor_type=3)
    comp1 = XOR.closeCells(comp1)  # close cells (adding point on lateral faces)
    #C.convertPyTree2File(comp1, "canal_adapted.cgns")
    return comp1

#==============================================================================
# periodicMeshAssembly

# GOAL : Conformal assembly of the input components

# IN: t1             : 3D NGON mesh, first operand
# IN: t2             : 3D NGON mesh, second operand
# IN: truezone1      : zone number in t1 which is the real component
# IN: truezone1      : zone number in t2 which is the real component

# OUT: returns the conformal assembly
#==============================================================================
def periodicMeshAssembly(t1, t2, TOL, real_zone_list):

    t0 = time.time()
    assembly = XOR.booleanUnion(t1, t2, tol=0., jtol=TOL, agg_mode=2, improve_qual=1, multi_zone=True)
    t1 = time.time()
    print (' - Boolean CPU time : ',t1-t0,'s')
    #C.convertPyTree2File(assembly, 'assembly.cgns')

    print("Extract and merge true components ...")
    ##### FUSION EN UN SEUL BLOC DES PARTIES UTILES
    zones = I.getZones(assembly)
    zs = []
    sz = len(real_zone_list)
    for i in range(sz):
        zs.append(zones[real_zone_list[i]])

    assembly = XOR.concatenate(zs, tol = TOL) # now one single block

    print("Check conformity ...")
    ##### VERIFICATION 1 : CONFORMITE (PAS DE FACES INTERNES)
    valid = XOR.isConformalNConnex(assembly, 1)
    if valid:
        print('Boolean OK : Conformal')
    else:
        print('Boolean ERROR : Internal non conformities ')
        ext = P.exteriorFaces(assembly)
        C.convertPyTree2File(ext, 'ext.cgns')
        import sys; sys.exit()

    return assembly


#==============================================================================
# detectMatchPerioAnomaly

# GOAL : Conformal assembly of the input components

# IN: t             : 3D NGON mesh
# IN: nperio        : nb of periodic boundaries (tipically 2 for a sector)
# IN: ncontours     : nb of closed loops in the periodic boundaries

# OUT: returns the conformal assembly
#==============================================================================
def detectMatchPerioAnomaly(t, nperio, ncontours):
    bcs = C.extractBCOfType(t, 'BCMatch')
    if bcs == []:
        print('No BC Match in mesh')
        return True

    bcs = T.join(bcs)
    bcs = T.splitConnexity(bcs)

    zs = I.getZones(bcs)
    if len(zs) != nperio :
        print ('BCMatchPerio ERROR : more than 2 zones')
        C.convertPyTree2File(zs, 'badperio.cgns')
        return True

    bcs = T.join(bcs)
    ext = P.exteriorFaces(bcs)
    #C.convertPyTree2File(ext, 'exf.cgns')
    ext = T.splitConnexity(ext)
    #C.convertPyTree2File(ext, 'split.cgns')
    zs = I.getZones(ext)
    if len(zs) != ncontours :
        print ('BCMatchPerio ERROR : missing some perio pairs')
        C.convertPyTree2File(zs, 'badperio.cgns')
        return True

    return False # still need to double check, not a complete checking..
