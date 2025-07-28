import XCore.PyTree as XC
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Transform.PyTree as T
import Generator.PyTree as G
import Post.PyTree as P
from Converter.Internal import E_NpyInt as E_NpyInt
import ToolboxIBM_CODA_MPI as TIBM
import numpy

TOL = 1e-9
def createQuadSurfaceFromNgonPointListSmallFace(a,PL):

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
            indices_face_5nodes = findFifthNode(idx,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a)
        for i in range(len_face):
            if len_face==4 or (len_face!=4 and not i in indices_face_5nodes):
                new_EC_ngon_faces[idx_new] = EC_faces[idx_face_init+i]
                idx_new += 1

    zone = Internal.newZone(name="Zone",zsize=[[nb_vertices_a,len_new_faces,0]],ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=zone)
    Internal.newDataArray('CoordinateX', value=coordsx_a, parent=gc)
    Internal.newDataArray('CoordinateY', value=coordsy_a, parent=gc)
    Internal.newDataArray('CoordinateZ', value=coordsz_a, parent=gc)
    elt = Internal.newElements(name="GridElements", etype=7, econnectivity=new_EC_ngon_faces, erange=[1, len_new_faces], eboundary=0, parent=zone)

    return zone

def createQuadSurfaceFromNgonPointListBigFace(a,cranges,indices_owners=[],dimPb=3):

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
                        indices_face_5nodes = findFifthNode(idx_face,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a)
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
    gcBI = Internal.newGridCoordinates(parent=zsnc_big_internal)
    Internal.newDataArray('CoordinateX', value=coordsx_a, parent=gcBI)
    Internal.newDataArray('CoordinateY', value=coordsy_a, parent=gcBI)
    Internal.newDataArray('CoordinateZ', value=coordsz_a, parent=gcBI)
    Internal.newElements(name="GridElements", etype=7, econnectivity=flat_list_bigfaces, erange=[1, len_new_faces], eboundary=0, parent=zsnc_big_internal)
    return zsnc_big_internal

def _createBCNearMatch(a_hexa,bcnearmatch):

    z = Internal.getZones(a_hexa)
    if Internal.getValue(bcnearmatch)[0][1]>0:
        f = P.exteriorFaces(z[0])
        hook = C.createHook(f,"elementCenters")
        ids = C.identifyElements(hook, bcnearmatch, tol=TOL)

        ids = ids[ids[:] > -1] #solo positivi
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))] #semplicemente ids-1
        if len(ids)>0:
            zf = T.subzone(f,ids, type='elements')
            TIBM._addBC2ZoneLoc(z[0], "QuadNQuad", "FamilySpecified:QuadNQuad", zf)
        C.freeHook(hook)
    return None

def _createBCStandard(a_hexa,a):
    nodes_bcs = Internal.getNodesFromType(a,"BC_t")
    for node_bc in nodes_bcs:
        bctype = Internal.getValue(node_bc)
        bcname = Internal.getName(node_bc)
        PL = Internal.getNodeFromName(node_bc,"PointList")[1]
        _createQuadConnectivityFromNgonPointList(a_hexa, a, PL, bcname, bctype)
    return None

def findFifthNode(idx,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a):
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
    point5_to_remove = numpy.where(numpy.abs(dists_from_center)<0.8*max_dist)[0]#[0]

    return point5_to_remove

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
        pos8 = numpy.where(EC==point8)[0][0]
        pos_corner_point = (pos8 + 2)%4
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

def createEmptyQuadZone():
    zone = Internal.newZone(name="empty",zsize=[[0,0]],ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=zone)
    Internal.newDataArray('CoordinateX', value=numpy.empty(0), parent=gc)
    Internal.newDataArray('CoordinateY', value=numpy.empty(0), parent=gc)
    Internal.newDataArray('CoordinateZ', value=numpy.empty(0), parent=gc)
    elt = Internal.newElements(name="GridElements", etype=7, econnectivity=numpy.empty(0), erange=[1, 0], eboundary=0, parent=zone)
    return zone


def _createQuadConnectivityFromNgonPointList(a_hexa, a, PL, bcname, bctype):

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
            indices_face_5nodes = findFifthNode(idx,offset_faces,length_faces,EC_faces,coordsx_a,coordsy_a,coordsz_a)
        for i in range(len_face):
            if len_face==4 or (len_face!=4 and not i in indices_face_5nodes):
                new_EC_ngon_faces[idx_new] = EC_faces[idx_face_init+i]
                idx_new += 1

    elts = Internal.getNodesFromType(a_hexa,"Elements_t")
    maxElt = Internal.getNodeFromName(elts[-1],"ElementRange")[1][1]
    zones = Internal.getZones(a_hexa)
    Internal.newElements(name=bcname, etype=7, econnectivity=new_EC_ngon_faces, erange=[maxElt+1, maxElt+len_new_faces], eboundary=0, parent=zones[0])
    C._addBC2Zone(zones[0],bcname,bctype, elementRange=[maxElt+1,maxElt+len_new_faces])
    return None
