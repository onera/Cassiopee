# IMPROVEMENTS
# cartRX in offsets ? 
# snear to adjust as in FAST/IBM
import Converter.PyTree as C
import Converter.Filter2 as Filter2
import Converter.Mpi as Cmpi
import Transform.PyTree as T 
import Converter.Internal as Internal
import Connector.PyTree as X
import Dist2Walls.PyTree as DTW
import Post.PyTree as P
import Geom.IBM as D_IBM
import XCore.PyTree as XC
from . import PyTree as G
import numpy
from Converter.Internal import E_NpyInt as E_NpyInt
__TOL__ = 1e-9

# Generation of the list of offset surfaces starting from tb
# IN: offsetValues : list of float values defining the offset distance to tb
# returns a tree toffset
def generateListOfOffsets__(tb, offsetValues=[], dim=3):
    if offsetValues==[]: return []

    if Cmpi.rank==0: print('Generating list of offsets...start',flush=True)

    BB = G.bbox(tb)
    ni = 150; nj = 150; nk = 150
    XRAYDIM1 = 3*ni; XRAYDIM2 = 3*nj

    offsetValMin = min(offsetValues)
    offsetValMax = max(offsetValues)

    alpha=1.1
    delta = alpha*offsetValMax
    xmin = BB[0]-delta; ymin = BB[1]-delta; zmin = BB[2]-delta
    xmax = BB[3]+delta; ymax = BB[4]+delta; zmax = BB[5]+delta
    # CARTRX
    bb_tb = G.bbox(tb)  
    xmin_core = bb_tb[0]-2*offsetValMin;  
    ymin_core = bb_tb[1]-2*offsetValMin;  
    zmin_core = bb_tb[2]-2*offsetValMin;  
    xmax_core = bb_tb[3]+2*offsetValMin;  
    ymax_core = bb_tb[4]+2*offsetValMin;  
    zmax_core = bb_tb[5]+2*offsetValMin;  

    hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1)
    if dim == 2:
        zmin = 0; zmax = 0
        zmin_core = 0.; zmax_core = 0.
        hk = 0.
    else:
        hk = (zmax-zmin)/(nk-1)    

    XC0 = (xmin_core, ymin_core, zmin_core); XF0 = (xmin, ymin, zmin)
    XC1 = (xmax_core, ymax_core, zmax_core); XF1 = (xmax, ymax, zmax)    
    b = G.cartRx3(XC0, XC1, (hi,hj,hk), XF0, XF1, (1.2,1.2,1.2), dim=dim, rank=Cmpi.rank, size=Cmpi.size)
    #
    DTW._distance2Walls(b, tb, type='mininterf', loc='nodes', signed=0)
    bodies = []
    for baseB in Internal.getBases(tb):
        bodies.append(Internal.getZones(baseB))
    nbodies = len(bodies)
    BM = numpy.ones((1, nbodies), dtype=numpy.int32)
    t = C.newPyTree(["BASE",Internal.getZones(b)])
    X._blankCells(t, bodies, BM, blankingType='node_in',dim=dim,XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2)
    C._initVars(t,'{TurbulentDistance}={TurbulentDistance}*({cellN}>0.)-{TurbulentDistance}*({cellN}<1.)')

    t_offset = C.newPyTree()
    for no_offset, offsetval in enumerate(offsetValues):
        #print(offsetval, flush=True)
        iso = P.isoSurfMC(t, 'TurbulentDistance',offsetval)
        iso = Cmpi.allgatherZones(iso)
        iso = C.convertArray2Tetra(iso)
        iso = T.join(iso)
        iso[0]='z_offset_%d'%(no_offset)
        C._addBase2PyTree(t_offset, 'OFFSET_%d'%(no_offset))
        t_offset[2][no_offset+1][2]=[iso]
    return t_offset

# Generates an isotropic skeleton mesh to be adapted then by AMR
def generateSkeletonMesh__(tb, dim=3, levelSkel=6):
    octreeMode = 0 # 1 : adjust of dfar
    surfaces=[]; dfarList=[]; snearsList=[]
    # list of dfars
    bodies = Internal.getZones(tb)
    for c, z in enumerate(bodies):
        n = Internal.getNodeFromName2(z, 'dfar')
        n2 = Internal.getNodeFromName2(z, 'snear')
        if n is not None: 
            dfarloc = Internal.getValue(n)
            snearloc = Internal.getValue(n2)
            if dfarloc > -1:#body snear is only considered if dfar_loc > -1
                surfaces.append(z)
                snearloc = 2**levelSkel*snearloc
                snearsList.append(snearloc)
                dfarList.append(dfarloc)

    o = G.octree(tb, snearList=snearsList, dfarList=dfarList, balancing=1, octreeMode=octreeMode)
    refined=True
    G._getVolumeMap(o)
    volminAll = C.getMinValue(o,"centers:vol")
    tol_vol = 1e-2*volminAll
    while refined:
        C._initVars(o,'{centers:indicator}=({centers:vol}>%g)'%(volminAll+tol_vol))
        if C.getMaxValue(o,"centers:indicator")==1.:
            o = G.adaptOctree(o, 'centers:indicator', balancing=1)
            G._getVolumeMap(o)
        else: 
            refined=False
            break
    C._rmVars(o, ['centers:indicator','centers:vol'])
   
    if dim==2: T._addkplane(o)
    o = C.convertArray2NGon(o)
    o = G.close(o)
    addPhysicalBCs__(o, tb, dim=dim)
    Internal._adaptNGon32NGon4(o)
    return o

def tagInsideOffset__(o, offset1=None, offset2=None, dim=3, h_target=-1.):
    to = C.newPyTree(["OCTREE"]); to[2][1][2]=Internal.getZones(o)
    C._initVars(to,'centers:indicator',0.)

    if dim==2:
        offset1 = T.addkplane(offset1)
        offset2 = T.addkplane(offset2)       
    bodies1 = [Internal.getZones(offset1)] # tag for refinement outside of offset1
    bodies2 = [Internal.getZones(offset2)] # tag for refinement inside of offset2 
    
    BM = numpy.ones((1,1),dtype=Internal.E_NpyInt)
    XRAYDIM1 = 400; XRAYDIM2 = 400

    C._initVars(to, "cellNOut",1.)
    C._initVars(to, "cellNIn",1.)   

    to = X.blankCells(to, bodies1, BM, blankingType='node_in', 
                      XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=dim, 
                      cellNName='cellNOut')
    to = X.setHoleInterpolatedPoints(to, depth=-1, cellNName='cellNOut', loc='nodes')
    to = X.blankCells(to, bodies2, BM, blankingType='node_in', 
                      XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=dim, 
                      cellNName='cellNIn')

    C._initVars(to,'{cellN}=({cellNIn}<1)*({cellNOut}>0.)')
    to = C.node2Center(to,["cellN"])
    C._initVars(to,"{centers:indicator}=({centers:cellN}>0.)")
    #
    G._getVolumeMap(to)
    vol_target = h_target**dim * 1.01 # plus a tolerance
    C._initVars(to,"{centers:indicator}={centers:indicator}*({centers:vol}>%g)"%vol_target)
    #
    C._rmVars(to, ["cellN","cellNIn","cellNOut","centers:cellN","centers:vol","centers:h"])
    o = Internal.getZones(to)[0]
    return o

def createQuadSurfaceFromNgonPointListBigFace__(a, cranges, indices_owners=[], dimPb=3):
    faces = Internal.getNodeFromName(a,"NGonElements")
    vol_cells = Internal.getNodeFromName(a,"NFaceElements")

    EC_volcells = Internal.getNodeFromName(vol_cells,"ElementConnectivity")[1]
    EC_faces = Internal.getNodeFromName(faces,"ElementConnectivity")[1]

    ER_volcells = Internal.getNodeFromName(vol_cells,"ElementRange")[1]
    ER_faces = Internal.getNodeFromName(faces,"ElementRange")[1]

    offset_volcells = Internal.getNodeFromName(vol_cells,"ElementStartOffset")[1]
    offset_faces = Internal.getNodeFromName(faces,"ElementStartOffset")[1]

    length_volcells = offset_volcells[1:] - offset_volcells[:-1]
    length_faces = offset_faces[1:] - offset_faces[:-1]
    lens_vol = numpy.unique(length_volcells)
    lens_faces = numpy.unique(length_faces)
    coordsx_a = Internal.getNodeFromName(a,"CoordinateX")[1]
    coordsy_a = Internal.getNodeFromName(a,"CoordinateY")[1]
    coordsz_a = Internal.getNodeFromName(a,"CoordinateZ")[1]
    nb_vertices_a = len(coordsx_a)

    if len(indices_owners) == 0:
        indices_owners = range(len(length_volcells))

    if dimPb == 3: n_smallfaces = 4
    else: n_smallfaces = 2

    list_bigface = []
    for index_volume in indices_owners:
        stride = cranges[index_volume]
        stride_offset = numpy.cumsum(stride)
        idx_vol_init = offset_volcells[index_volume]

        idx_sides = []
        for i in range(6):
            if stride[i] == n_smallfaces:
                idx_sides.append(i)
        if idx_sides != []:
            for idx_side in idx_sides:
                indices_faces = []
                start_idx = stride_offset[idx_side-1] if idx_side>0 else 0
                for i in range(start_idx,stride_offset[idx_side]):
                    face = EC_volcells[idx_vol_init+i]
                    indices_faces.append(face)
                conn_Nfaces = numpy.zeros((n_smallfaces,4),dtype=E_NpyInt)
                for i,idx_face in enumerate(indices_faces):
                    idx_face_init = offset_faces[idx_face-1]
                    len_face = length_faces[idx_face-1]

                    if len_face != 4:
                        indices_face_5nodes = findFifthNode__(idx_face,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a)
                    k = 0
                    for j in range(len_face):
                        if len_face==4 or (len_face!=4 and not j in indices_face_5nodes):
                            conn_Nfaces[i][k] = EC_faces[idx_face_init+j]
                            k+=1
                if dimPb == 3: [point0,point1,point2,point3] = reorderNodesInCanonicalOrderForBigFace3D(conn_Nfaces)
                else: [point0,point1,point2,point3] = reorderNodesInCanonicalOrderForBigFace2D(conn_Nfaces)
                list_bigface.append([point0,point1,point2,point3])

    len_new_faces = len(list_bigface)
    flat_list_bigfaces = numpy.array([x for xs in list_bigface for x in xs])

    zsnc_big_internal = Internal.newZone(name="Zone",zsize=[[nb_vertices_a,len_new_faces,0]],ztype="Unstructured")
    gcBI = Internal.newGridCoordinates(parent = zsnc_big_internal)
    Internal.newDataArray('CoordinateX', value = coordsx_a, parent = gcBI)
    Internal.newDataArray('CoordinateY', value = coordsy_a, parent = gcBI)
    Internal.newDataArray('CoordinateZ', value = coordsz_a, parent = gcBI)
    Internal.newElements(name = "GridElements", etype = 7, econnectivity = flat_list_bigfaces, erange = [1, len_new_faces], eboundary = 0, parent = zsnc_big_internal)
    return zsnc_big_internal

def createQuadSurfaceFromNgonPointListSmallFace__(a,PL):
    faces = Internal.getNodeFromName(a,"NGonElements")
    vol_cells = Internal.getNodeFromName(a,"NFaceElements")

    EC_volcells = Internal.getNodeFromName(vol_cells,"ElementConnectivity")[1]
    EC_faces = Internal.getNodeFromName(faces,"ElementConnectivity")[1]

    ER_volcells = Internal.getNodeFromName(vol_cells,"ElementRange")[1]
    ER_faces = Internal.getNodeFromName(faces,"ElementRange")[1]

    offset_volcells = Internal.getNodeFromName(vol_cells,"ElementStartOffset")[1]
    offset_faces = Internal.getNodeFromName(faces,"ElementStartOffset")[1]

    length_volcells = offset_volcells[1:] - offset_volcells[:-1]
    length_faces = offset_faces[1:] - offset_faces[:-1]
    lens_vol = numpy.unique(length_volcells)
    lens_faces = numpy.unique(length_faces)

    coordsx_a = Internal.getNodeFromName(a,"CoordinateX")[1]
    coordsy_a = Internal.getNodeFromName(a,"CoordinateY")[1]
    coordsz_a = Internal.getNodeFromName(a,"CoordinateZ")[1]
    nb_vertices_a = len(coordsx_a)

    len_new_faces = len(PL)
    len_EC_faces = len(PL)*4
    new_EC_ngon_faces = numpy.zeros(len_EC_faces,dtype=E_NpyInt)
    idx_new = 0
    for idx in PL :
        idx_face_init = offset_faces[idx-1]
        len_face = length_faces[idx-1]
        if len_face !=4:
            indices_face_5nodes = findFifthNode__(idx,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a)
        for i in range(len_face):
            if len_face==4 or (len_face!=4 and not i in indices_face_5nodes):
                new_EC_ngon_faces[idx_new] = EC_faces[idx_face_init+i]
                idx_new += 1

    zone = Internal.newZone(name="Zone",zsize=[[nb_vertices_a,len_new_faces,0]],ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent = zone)
    Internal.newDataArray('CoordinateX', value = coordsx_a, parent = gc)
    Internal.newDataArray('CoordinateY', value = coordsy_a, parent = gc)
    Internal.newDataArray('CoordinateZ', value = coordsz_a, parent = gc)
    elt = Internal.newElements(name = "GridElements", etype = 7, econnectivity = new_EC_ngon_faces, erange = [1, len_new_faces], eboundary = 0, parent = zone)
    return zone

def findFifthNode__(idx,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a):
    idx_face_init = offset_faces[idx-1]
    len_face = length_faces[idx-1]
    EC_ngon_faces = EC_faces[idx_face_init:idx_face_init+len_face]
    x_5 = coordsx_a[EC_ngon_faces-1]
    y_5 = coordsy_a[EC_ngon_faces-1]
    z_5 = coordsz_a[EC_ngon_faces-1]

    x_center = (numpy.max(x_5)+numpy.min(x_5))/2
    y_center = (numpy.max(y_5)+numpy.min(y_5))/2
    z_center = (numpy.max(z_5)+numpy.min(z_5))/2

    dists_from_center = ((x_5-x_center)**2+(y_5-y_center)**2+(z_5-z_center)**2)**0.5
    max_dist = numpy.max(dists_from_center)
    point5_to_remove = numpy.where(numpy.abs(dists_from_center)<0.8*max_dist)[0]
    return point5_to_remove
 
def createPseudoBCQuadNQuadIntra__(a, owners, levels, halo_levels, neighbours, cranges, dimPb=3):
    point_list_small = []
    point_list_big = []
    indices_owners = []

    for PL in range(len(owners)):
        own = owners[PL]
        nei = neighbours[PL]
        own_lev = levels[own]
        nei_lev = levels[nei]
        if nei > -1:
            if own_lev != nei_lev:
                point_list_small.append(PL+1)
                point_list_big.append(PL+1)
                indices_owners.append(own)
                indices_owners.append(nei)
    if len(point_list_small) > 0:
        zsnc_small_internal = createQuadSurfaceFromNgonPointListSmallFace__(a, point_list_small)
    else:
        zsnc_small_internal = createEmptyQuadZone__()
    if len(point_list_big) > 0:
        zsnc_big_internal = createQuadSurfaceFromNgonPointListBigFace__(a, cranges, indices_owners, dimPb=dimPb)
    else:
        zsnc_big_internal = createEmptyQuadZone__()
    zone_nonconformal_in_proc = T.join(zsnc_big_internal, zsnc_small_internal)
    return zone_nonconformal_in_proc

def createPseudoBCQuadNQuadInter__(a, owners, levels, halo_levels, neighbours, cranges, dimPb=3):
    point_list_commpatch = []
    zgc_nodes = Internal.getNodesFromType(a,"GridConnectivity1to1_t")
    for zgc_node in zgc_nodes:
        PL = Internal.getNodeFromName(zgc_node,"PointList")[1]
        point_list_commpatch.append(PL)

    point_list_small = []
    point_list_big = []
    for nei_proc in range(len(zgc_nodes)):
        PL_nei_proc = point_list_commpatch[nei_proc]
        for i in range(len(PL_nei_proc)):
            PL = PL_nei_proc[i]
            own = owners[PL-1]
            own_lev = levels[own]
            nei_lev = halo_levels[nei_proc][i]
            if own_lev < nei_lev:
                point_list_big.append(PL)
            elif own_lev > nei_lev:
                point_list_small.append(PL)
    if len(point_list_small)>0:
        zsnc_small_MPI = createQuadSurfaceFromNgonPointListSmallFace__(a, point_list_small)
    else:
        zsnc_small_MPI = createEmptyQuadZone__()

    if len(point_list_big)>0:
        point_list_big = numpy.array(point_list_big)
        indices_owners = owners[point_list_big-1]
        zsnc_big_MPI = createQuadSurfaceFromNgonPointListBigFace__(a, cranges, indices_owners, dimPb)
    else:
        zsnc_big_MPI = createEmptyQuadZone__()
  
    zone_nonconformal_between_procs = T.join(zsnc_big_MPI, zsnc_small_MPI)    
    return zone_nonconformal_between_procs

def reorderNodesInCanonicalOrderForBigFace2D(conn_2faces):
    point0=None; point1=None;point2=None;point3=None
    unique, counts = numpy.unique(conn_2faces, return_counts=True)
    indices_count_2 = numpy.argwhere(counts==2)
    indices_count_2 = indices_count_2.reshape(2)
    indices_count_1 = numpy.argwhere(counts==1)
    indices_count_1 = indices_count_1.reshape(4)
    points45 = unique[indices_count_2]
    corners_sort = unique[indices_count_1]
    point0 = corners_sort[0] # small index
    if point0 in conn_2faces[0]:
        face0 = 0
        face1 = 1
    else:
        face0 = 1
        face1 = 0
    conn_face0 = conn_2faces[face0]
    pos_point0 = numpy.where(conn_face0==point0)[0][0]

    if conn_face0[(pos_point0+1)%4] in points45:
        point4 = conn_face0[(pos_point0+1)%4]
        point5 = conn_face0[(pos_point0+2)%4]
        point3 = conn_face0[(pos_point0+3)%4]
    else:
        point4 = conn_face0[(pos_point0-1)%4]
        point5 = conn_face0[(pos_point0-2)%4]
        point3 = conn_face0[(pos_point0-3)%4]

    conn_face1 = conn_2faces[face1]
    pos_point4 = numpy.where(conn_face1==point4)[0][0]

    if conn_face1[(pos_point4+1)%4] in points45:
        point2 = conn_face1[(pos_point4+2)%4]
        point1 = conn_face1[(pos_point4+3)%4]
    else:
        point2 = conn_face1[(pos_point4-2)%4]
        point1 = conn_face1[(pos_point4+1)%4]
        
    if point1>point3:
        tmp = point1
        point1 = point3
        point3 = tmp
    return [point0,point1,point2,point3]

def reorderNodesInCanonicalOrderForBigFace3D(conn_4faces):
    point0=None; point1=None;point2=None;point3=None
    unique, counts = numpy.unique(conn_4faces, return_counts=True)
    index_count_4 = numpy.argwhere(counts==4)[0][0]
    indices_count_1 = numpy.argwhere(counts==1)
    indices_count_1 = indices_count_1.reshape(4)
    point8 = unique[index_count_4]
    corners = unique[indices_count_1]

    corner_points = []
    pos_corner_points = []
    for EC in conn_4faces:
        pos8 = numpy.where(EC == point8)[0][0]
        pos_corner_point = (pos8 + 2) % 4
        pos_corner_points.append(pos_corner_point)
        corner_points.append(EC[pos_corner_point])

    elt0 = numpy.argmin(corner_points)
    point0 = corner_points[elt0]
    pos_corner_in_elt0 = pos_corner_points[elt0]

    point4 = conn_4faces[elt0][(pos_corner_in_elt0+1)%4]
    point7 = conn_4faces[elt0][(pos_corner_in_elt0-1)]

    elt1 = numpy.where(conn_4faces==point4)[0]
    elt1 = elt1[elt1!=elt0][0]

    point1 = corner_points[elt1]
    pos_corner_in_elt1 = pos_corner_points[elt1]

    point5 = list(set(conn_4faces[elt1])-set([point4,point8,point1]))[0]

    elt2 = numpy.where(conn_4faces==point5)[0]
    elt2 = elt2[elt2!=elt1][0]

    point2 = corner_points[elt2]
    pos_corner_in_elt2 = pos_corner_points[elt2]

    point6 = list(set(conn_4faces[elt2])-set([point5,point8,point2]))[0]
    elt3 = list(set([0,1,2,3])-set([elt0,elt1,elt2]))[0]
    point3 = list(set(conn_4faces[elt3])-set([point6,point7,point8]))[0]
    return [point0,point1,point2,point3]

def createEmptyQuadZone__():
    zone = Internal.newZone(name="empty",zsize=[[0,0]],ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent = zone)
    Internal.newDataArray('CoordinateX', value = numpy.empty(0), parent = gc)
    Internal.newDataArray('CoordinateY', value = numpy.empty(0), parent = gc)
    Internal.newDataArray('CoordinateZ', value = numpy.empty(0), parent = gc)
    elt = Internal.newElements(name = "GridElements", etype = 7, econnectivity = numpy.empty(0), erange = [1, 0], eboundary = 0, parent = zone)
    return zone

def _createQuadConnectivityFromNgonPointList__(a_hexa, a, PL, bcname, bctype):
    faces = Internal.getNodeFromName(a, "NGonElements")
    vol_cells = Internal.getNodeFromName(a, "NFaceElements")

    EC_volcells = Internal.getNodeFromName(vol_cells, "ElementConnectivity")[1]
    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]

    ER_volcells = Internal.getNodeFromName(vol_cells, "ElementRange")[1]
    ER_faces = Internal.getNodeFromName(faces, "ElementRange")[1]

    offset_volcells = Internal.getNodeFromName(vol_cells, "ElementStartOffset")[1]
    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]

    length_volcells = offset_volcells[1:] - offset_volcells[:-1]
    length_faces = offset_faces[1:] - offset_faces[:-1]
    lens_vol = numpy.unique(length_volcells)
    lens_faces = numpy.unique(length_faces)

    coordsx_a = Internal.getNodeFromName(a, "CoordinateX")[1]
    coordsy_a = Internal.getNodeFromName(a, "CoordinateY")[1]
    coordsz_a = Internal.getNodeFromName(a, "CoordinateZ")[1]
    nb_vertices_a = len(coordsx_a)

    len_new_faces = len(PL)
    len_EC_faces = len(PL) * 4
    new_EC_ngon_faces = numpy.zeros(len_EC_faces, dtype=E_NpyInt)
    idx_new = 0
    for idx in PL:
        idx_face_init = offset_faces[idx - 1]
        len_face = length_faces[idx - 1]
        if len_face != 4:
            indices_face_5nodes = findFifthNode__(idx, offset_faces, length_faces, EC_faces, coordsx_a, coordsy_a, coordsz_a)
        for i in range(len_face):
            if len_face == 4 or (len_face != 4 and not i in indices_face_5nodes):
                new_EC_ngon_faces[idx_new] = EC_faces[idx_face_init + i]
                idx_new += 1

    elts = Internal.getNodesFromType(a_hexa, "Elements_t")
    maxElt = Internal.getNodeFromName(elts[-1], "ElementRange")[1][1]
    zones = Internal.getZones(a_hexa)
    Internal.newElements(
        name=bcname,
        etype=7,
        econnectivity=new_EC_ngon_faces,
        erange=[maxElt + 1, maxElt + len_new_faces],
        eboundary=0,
        parent=zones[0]
    )
    C._addBC2Zone(zones[0], bcname, bctype, elementRange=[maxElt + 1, maxElt + len_new_faces])
    return None
 
def _createBCNearMatch__(a_hexa, bcnearmatch):
    z = Internal.getZones(a_hexa)
    if Internal.getValue(bcnearmatch)[0][1]>0:
        f = P.exteriorFaces(z[0])
        hook = C.createHook(f,"elementCenters")
        ids = C.identifyElements(hook, bcnearmatch, tol=__TOL__)

        ids = ids[ids[:] > -1] 
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))] 
        if len(ids)>0:
            zf = T.subzone(f,ids, type='elements')
            _addBC2Zone__(z[0], "QuadNQuad", "FamilySpecified:QuadNQuad", zf)
        C.freeHook(hook)
    return None

def _createBCStandard__(a_hexa, a):
    nodes_bcs = Internal.getNodesFromType(a,"BC_t")
    for node_bc in nodes_bcs:
        bctype = Internal.getValue(node_bc)
        bcname = Internal.getName(node_bc)
        PL = Internal.getNodeFromName(node_bc,"PointList")[1]
        _createQuadConnectivityFromNgonPointList__(a_hexa, a, PL, bcname, bctype)
    return None

def adaptMesh__(FILE_SKELETON, hmin, tb, bbo, toffset=None, dim=3):
    o, res = XC.loadAndSplitNGon(FILE_SKELETON)
    Cmpi.barrier()
    gcells = res[5]
    gfaces = res[6]
    comm = res[1]
    if dim==3: normal2D=None
    else: normal2D = numpy.array([0.0,0.0,1.0])
    hookAM = XC.AdaptMesh_Init(o, normal2D, comm=comm, gcells=gcells, gfaces=gfaces)

    offset_zones = Internal.getZones(toffset)
    noffsets = len(offset_zones)
    offset_inside = Internal.getZones(tb)
    for i in range(noffsets-1, -1,-1):        
        #XC.AdaptMesh_LoadBalance(hookAM)       
        offsetloc = offset_zones[i]
        hx = hmin * 2**i 

        adapting=True
        while adapting:        
            o = tagInsideOffset__(o, offset1=offset_inside, offset2=offsetloc, dim=dim, h_target=hx)
            indicMax = C.getMaxValue(o,"centers:indicator")
            indicMax = Cmpi.allgather(indicMax)
            indicMax = max(indicMax)
            if indicMax<1.:
                adapting=False
                C._rmVars(o,["centers:indicator"])
                break
            else:
                f = Internal.getNodeFromName(o, 'indicator')[1]
                REF = f.astype(dtype=Internal.E_NpyInt)
                XC.AdaptMesh_AssignRefData(hookAM, REF)
                #XC.AdaptMesh_LoadBalance(hookAM)       
                XC.AdaptMesh_Adapt(hookAM)
                o = XC.AdaptMesh_ExtractMesh(hookAM, conformize=1)
                o = Internal.getZones(o)[0]

    o = XC.AdaptMesh_ExtractMesh(hookAM, conformize=1)
    owners = XC.AdaptMesh_ExtractOwners(hookAM)
    levels = XC.AdaptMesh_ExtractCellLevels(hookAM)
    halo_levels = XC.AdaptMesh_ExtractHaloCellLevels(hookAM)
    neighbours = XC.AdaptMesh_ExtractNeighbours(hookAM)
    cranges = XC.AdaptMesh_ExtractCellRanges(hookAM)    
    
    cart_hexa = XC.AdaptMesh_ExtractMesh(hookAM, conformize=0)
    XC.AdaptMesh_Exit(hookAM)
    cart_hexa = Internal.getZones(cart_hexa)[0] 
        
    zone_nonconformal_inter = createPseudoBCQuadNQuadInter__(o, owners, levels, halo_levels, neighbours, cranges, dimPb=dim)
    zone_nonconformal_intra = createPseudoBCQuadNQuadIntra__(o, owners, levels, halo_levels, neighbours, cranges, dimPb=dim)
    zone_nonconformal = T.join(zone_nonconformal_inter, zone_nonconformal_intra) 

    _createBCNearMatch__(cart_hexa,zone_nonconformal)
    _createBCStandard__(cart_hexa, o)
    del o

    Cmpi._setProc(cart_hexa, Cmpi.rank)  
    return cart_hexa

def _addBC2Zone__(z, bndName, bndType, zbc, loc='FaceCenter', zdnrName=None):
    s = bndType.split(':')
    bndType1 = s[0]
    if len(s) > 1: bndType2 = s[1]
    else: bndType2 = ''

    # Analyse zone zbc
    dims = Internal.getZoneDim(zbc)
    neb = dims[2] # nbre d'elts de zbc

    eltType, nf = Internal.eltName2EltNo(dims[3]) # type d'elements de zbc
    # On cherche l'element max dans les connectivites de z
    maxElt = 0
    connects = Internal.getNodesFromType(z, 'Elements_t')
    for cn in connects:
        r = Internal.getNodeFromName1(cn, 'ElementRange')
        m = r[1][1]
        maxElt = max(maxElt, m)
    # on cree un nouveau noeud connectivite dans z1 (avec le nom de la zone z2)
    nebb = neb
    node = Internal.createUniqueChild(z, bndName, 'Elements_t', value=[eltType,nebb])
    Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                               value=[maxElt+1,maxElt+neb])
    oldc = Internal.getNodeFromName2(zbc, 'ElementConnectivity')[1]
    newc = numpy.copy(oldc)
    hook = C.createHook(z, 'nodes')
    ids = C.identifyNodes(hook, zbc)
    newc[:] = ids[oldc[:]-1]
    Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)

    zoneBC = Internal.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
    if len(s)==1:
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType)
    else: # familyspecified
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        Internal.createUniqueChild(info, 'FamilyName', 'FamilyName_t',
                                   value=bndType2)

    Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                                value='FaceCenter')
    Internal.createUniqueChild(info, 'ElementRange', 'IndexRange_t',
                                value=numpy.array([[maxElt+1,maxElt+neb]]))
    return None

#==================================================================
# Addition of Physical BCs on the skeleton NGON Mesh
#==================================================================
def addPhysicalBCs__(z_ngon, tb, dim=3):
    bbo = G.bbox(z_ngon)
    
    # determine where the symmetry plane is 
    baseSYM = Internal.getNodeFromName1(tb,"SYM")
    dir_sym = 0
    if baseSYM is not None:
        symplane = []
        for zsym in Internal.getZones(baseSYM):
            if C.getMaxValue(zsym,'centers:cellN')>0.:
                symplane.append(zsym)
        [xmin,ymin,zmin,xmax,ymax,zmax] = G.bbox(symplane)
        if abs(xmax-xmin) < __TOL__:
            if abs(xmin)< __TOL__: dir_sym=-1
            else: dir_sym=1
        elif abs(ymax-ymin) < __TOL__:
            if abs(ymin)< __TOL__: dir_sym=-2
            else: dir_sym=2            
        elif abs(zmax-zmin)<__TOL__:
            if abs(zmin)< __TOL__: dir_sym=-3
            else: dir_sym=3

    xmin = bbo[0]; ymin = bbo[1]; zmin = bbo[2]
    xmax = bbo[3]; ymax = bbo[4]; zmax = bbo[5]

    extFaces = P.exteriorFaces(z_ngon)
    extFaces = T.splitSharpEdges(extFaces, 80.)

    zbc_xmin=[]; zbc_xmax=[]; zbc_ymin=[]; zbc_ymax=[]; zbc_zmin=[]; zbc_zmax=[]
    for e in extFaces:
        bbe = G.bbox(e)
        xmine = bbe[0]; ymine = bbe[1]; zmine = bbe[2]
        xmaxe = bbe[3]; ymaxe = bbe[4]; zmaxe = bbe[5]
        dx = xmaxe-xmine; dy = ymaxe-ymine; dz = zmaxe-zmine
        if dx < __TOL__ and abs(xmine-xmin)<__TOL__: zbc_xmin.append(e)
        elif dx < __TOL__ and abs(xmaxe-xmax)<__TOL__: zbc_xmax.append(e)
        elif dy < __TOL__ and abs(ymine-ymin)<__TOL__: zbc_ymin.append(e)
        elif dy < __TOL__ and abs(ymaxe-ymax)<__TOL__: zbc_ymax.append(e)
        elif dz < __TOL__ and abs(zmine-zmin)<__TOL__: zbc_zmin.append(e)
        elif dz < __TOL__ and abs(zmaxe-zmax)<__TOL__: zbc_zmax.append(e)
    
    zbc_xmin = T.join(zbc_xmin); zbc_xmin[0]='xmin'
    zbc_xmax = T.join(zbc_xmax); zbc_xmax[0]='xmax'
    zbc_ymin = T.join(zbc_ymin); zbc_ymin[0]='ymin'
    zbc_ymax = T.join(zbc_ymax); zbc_ymax[0]='ymax'
    zbc_zmin = T.join(zbc_zmin); zbc_zmin[0]='zmin'
    zbc_zmax = T.join(zbc_zmax); zbc_zmax[0]='zmax'
    if dir_sym==-1:
        zbc_xmin[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane", subzone=zbc_xmin)
    else:        
        zbc_xmin[0] = C.getZoneName("BCFarfield")
        C._addBC2Zone(z_ngon, 'BCFarfield',"BCFarfield",subzone=zbc_xmin)

    if dir_sym== 1:
        zbc_xmax[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_xmax)
    else:        
        zbc_xmax[0] = C.getZoneName("BCFarfield")
        C._addBC2Zone(z_ngon, 'BCFarfield',"BCFarfield",subzone=zbc_xmax)
                       
    if dir_sym==-2:
        zbc_ymin[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_ymin)
    else:        
        zbc_ymin[0] = C.getZoneName("BCFarfield")
        C._addBC2Zone(z_ngon, 'BCFarfield',"BCFarfield",subzone=zbc_ymin)

    if dir_sym==2:
        zbc_ymax[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_ymax)
    else:        
        zbc_ymax[0] = C.getZoneName("BCFarfield")
        C._addBC2Zone(z_ngon, 'BCFarfield',"BCFarfield",subzone=zbc_ymax)
                
    if dim==2:
        zbc_zmin[0] = C.getZoneName("BCSymmetryPlane")
        zbc_zmax[0] = C.getZoneName("BCSymmetryPlane")              
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_zmin)
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_zmax)
    else:
        if dir_sym==-3:
            zbc_zmin[0] = C.getZoneName("BCSymmetryPlane")
            C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_zmin)
        else:        
            zbc_zmin[0] = C.getZoneName("BCFarfield")
            C._addBC2Zone(z_ngon, 'BCFarfield',"BCFarfield",subzone=zbc_zmin)

        if dir_sym==3:
            zbc_zmax[0] = C.getZoneName("BCSymmetryPlane")
            C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane",subzone=zbc_zmax)
        else:        
            zbc_zmax[0] = C.getZoneName("BCFarfield")
            C._addBC2Zone(z_ngon, 'BCFarfield',"BCFarfield",subzone=zbc_zmax)
    return None

#==================================================================
# Generation of the AMR mesh for IBMs 
#==================================================================
def generateAMRMesh(tb, levelMax=7, vmins=11, dim=3, check=False):
    FILE_SKELETON = 'skeleton.cgns'
    # levelSkel: initial refinement level of the skeleton octree
    # might be tuned 
    o = generateSkeletonMesh__(tb, dim=dim, levelSkel=levelMax)
    G._getVolumeMap(o)
    hmin_skel = (C.getMinValue(o,"centers:vol"))**(1/dim)
    hmin = hmin_skel * 2 ** (-levelMax)
    # mandatory save file for loadAndSplit for adaptation
    if Cmpi.rank==0:
        C.convertPyTree2File(o,FILE_SKELETON)
    Cmpi.barrier()  

    # layers of offsets that are inside the bbox of the domain (dfar_max)
    dfar_max = 0.
    for z in Internal.getZones(tb):
        dfarloc = Internal.getValue(Internal.getNodeFromName(z,"dfar"))
        dfar_max = max(dfarloc,dfar_max)

    bbo = G.bbox(o)
    dfarmax = min(bbo[3]-bbo[0], bbo[4]-bbo[1])
    if dim==3: dfarmax = min(dfarmax,bbo[5]-bbo[2])

    offsetValues = []
    offsetprev = 0.
    print(" Minimum spacing = ", hmin, flush=True)
    for no_adapt in range(len(vmins)):
        offsetloc = offsetprev+hmin*(2**no_adapt)*vmins[no_adapt]
        if offsetloc < 0.99*dfarmax:
            offsetValues.append(offsetloc)
            offsetprev=offsetloc    
    #generate list of offsets 
    print("Generate list of offsets for rank ", Cmpi.rank, flush=True)
    toff = generateListOfOffsets__(tb, offsetValues=offsetValues, dim=dim)
    if check and Cmpi.rank==0: C.convertPyTree2File(toff,"offset.cgns")

    # adaptation of the mesh wrt to the bodies (finest level) and offsets
    # only a part is returned per processor
    o = adaptMesh__(FILE_SKELETON, hmin, tb, bbo, toffset=toff, dim=dim)
    Cmpi._setProc(o, Cmpi.rank)
    t = C.newPyTree(['AMR',o])
    return t