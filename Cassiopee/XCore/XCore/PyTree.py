import Converter.Filter2 as Filter2
import Converter.Mpi as Cmpi
import Converter.Internal as I
import Converter.PyTree as C
from . import xcore

BCType_l = set(I.KNOWNBCS)

# IN : a : zone
# IN : cid : cell id as integer
def extractCell(a, cid):
    z = I.getZones(a)[0]
    m = C.getFields(I.__GridCoordinates__, z, api=3)[0]
    zmo = I.createZoneNode("cell"+str(cid), xcore.extractCell(m, cid))
    return zmo

# -- AdaptMesh_Init
# Initialize a structure for mesh adaptation
# IN: t : CGNS tree/base/zones, only a single zone is taken into account
# IN: normal2D: vecteur 2D
# IN: comm: tableaux de connectivités en parallèle issus de chunk2part
# IN: gcells: indices globaux des cellules en parallèle
# IN: gfaces: indices globaux des faces en parallèle
# OUT: hook on the structure required during the mesh adaptation process
def AdaptMesh_Init(t, normal2D=None, comm=[], gcells=None, gfaces=None):
    zones = I.getZones(t)
    if len(zones) != 1:
        print("WARNING: AdaptMesh_Init: only the first zone is taken into account.", flush=True)

    z = zones[0]
    array = C.getFields(I.__GridCoordinates__, z, api=3)[0]

    bcs = []
    zonebc = I.getNodeFromType1(z, 'ZoneBC_t')
    if zonebc is not None:
        zbc = I.getNodesFromType1(zonebc, 'BC_t')
        bc_count = 0
        for bc in zbc:
            plist = I.getNodeFromName1(bc, 'PointList')
            name = bc[0]
            #tag = I.getNodeFromName(bc, 'Tag')[1][0]
            try: tag = I.getNodeFromName1(bc, 'Tag')[1][0]
            except: tag = bc_count; bc_count += 1
            bctype = I.getValue(bc)
            bcs.append([plist[1], tag, name, bctype])

    return xcore.AdaptMesh_Init(array, normal2D, bcs, comm, gcells, gfaces)

# -- AdaptMesh_Exit
# free the memory allocated by AdaptMesh_Init
# IN: AM: hook created by AdaptMesh_Init
def AdaptMesh_Exit(AM):
    return xcore.AdaptMesh_Exit(AM)

# IN/OUT : AM: hook on the structure created by AdaptMesh_Init
# IN : REF: numpy of 0/1 integers used as the adaptation indicator for each leaf of the current adaptation data structure
# return None
def AdaptMesh_AssignRefData(AM, REF):
    return xcore.AdaptMesh_AssignRefData(AM, REF)

# load balance using ptsotch
# IN/OUT: AM: hook on the structure created by AdaptMesh_Init
def AdaptMesh_LoadBalance(AM):
    return xcore.AdaptMesh_LoadBalance(AM)

# adaptation of the mesh
# IN/OUT: hook on the structure refering to the adapted mesh
def AdaptMesh_Adapt(AM):
    return xcore.AdaptMesh_Adapt(AM)

# convert the hook into a mesh
# IN: AM: data structure
# IN: conformize: 0=HEXA, 1=NGON
# OUT: return a zone per proc with 1to1 grid connectivities and bcs
def AdaptMesh_ExtractMesh(AM, conformize=1, recoverBC=True):
    mesh, bcs, comm, procs = xcore.AdaptMesh_ExtractMesh(AM, conformize)
    name = 'Proc' + '%d'%Cmpi.rank
    zone = I.createZoneNode(name, mesh)

    if recoverBC is False:
        t = C.newPyTree([name, zone])
        return t

    if procs is not None and len(procs) > 0:
        I.newUserDefinedData(name='NeighbourProcessors', value=procs, parent=zone)

    if comm is not None:
        for data in comm:
            cont = I.createUniqueChild(zone, 'ZoneGridConnectivity', 'ZoneGridConnectivity_t')
            Name = 'Match_'+str(data[0])
            I.newGridConnectivity1to1(name=Name, donorName=str(data[0]), pointList=data[1], parent=cont)

    # add BCs
    if bcs is not None:
        for i in range(len(bcs)):
            bc = bcs[i]
            if len(bc) != 0:
                cont = I.createUniqueChild(zone, 'ZoneBC', 'ZoneBC_t')
                ptlist = bc[0]
                tag = bc[1]
                bcname = bc[2]
                bctype = bc[3]

                if bctype not in BCType_l:
                    bc = I.newBC(name=bcname, pointList=ptlist, family=bctype, parent=cont)
                else:
                    bc = I.newBC(name=bcname, pointList=ptlist, btype=bctype, parent=cont)

                I.newUserDefinedData(name='Tag', value=tag, parent=bc)

    t = C.newPyTree([name, zone])
    return t

def AdaptMesh_ExtractOwners(AM):
    return xcore.AdaptMesh_ExtractOwners(AM)

def AdaptMesh_ExtractNeighbours(AM):
    return xcore.AdaptMesh_ExtractNeighbours(AM)

def AdaptMesh_ExtractCellLevels(AM):
    return xcore.AdaptMesh_ExtractCellLevels(AM)

def AdaptMesh_ExtractCellRanges(AM):
    return xcore.AdaptMesh_ExtractCellRanges(AM)

def AdaptMesh_ExtractHaloCellLevels(AM):
    return xcore.AdaptMesh_ExtractHaloCellLevels(AM)


################################################################################

# Returns for each zone, exchanged fields
def exchangeFields(t, fldNames):
    zones = I.getZones(t)
    rfields = []
    for zone in zones:
        arr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
        pe = I.getNodeFromName(zone, 'ParentElements')
        if pe is None: raise ValueError('ParentElements not found.')
        fsolc = I.getNodeFromName2(zone, I.__FlowSolutionCenters__)
        if fsolc is None: raise ValueError('FlowSolutionCenters not found.')
        flds = []
        for fldName in fldNames:
            fld = I.getNodeFromName2(fsolc, fldName)
            if fld is None: raise ValueError(fldName, 'not found.')
            flds.append(fld[1])
        zgc = I.getNodeFromType1(zone, 'ZoneGridConnectivity_t')
        if zgc is None: raise ValueError('ZoneGridConnectivity not found')
        comms = I.getNodesFromType1(zgc, 'GridConnectivity1to1_t')
        if comms is None: raise ValueError('GridConnectivity1to1 not found')
        comm_list = []
        for comm in comms:
            nei_proc = int(I.getValue(comm))
            ptlist = I.getNodeFromName(comm, 'PointList')[1]
            comm_list.append([nei_proc, ptlist])
        rfields.append(xcore.exchangeFields(arr, pe[1], flds, comm_list))
    return rfields

def loadAndSplitElt(fileName):
    dt = Filter2.loadAsChunks(fileName)
    zones = I.getZones(dt)

    if len(zones) > 1:
        raise TypeError("loadAndSplitElt: one zone only.")

    z = zones[0]

    cx = I.getNodeFromName2(z, 'CoordinateX')[1]
    cy = I.getNodeFromName2(z, 'CoordinateY')[1]
    cz = I.getNodeFromName2(z, 'CoordinateZ')[1]

    XYZ = []
    XYZ.append(cx); XYZ.append(cy); XYZ.append(cz)

    cns = I.getNodesFromType1(z, 'Elements_t')

    chunks = []

    for cn in cns:
        name, stride = I.eltNo2EltName(cn[1][0])
        arr = I.getNodeFromName1(cn, 'ElementConnectivity')[1]
        chunks.append([name, stride, arr])

    parts = xcore.chunk2partElt(XYZ, chunks)

    zones = []

    for i, p in enumerate(parts):
        z = I.createZoneNode('Zone' + '%d'%Cmpi.rank + '_%d'%i, p)
        zones.append(z)

    t = C.newPyTree(['Base', zones])

    Cmpi._setProc(t, Cmpi.rank)

    return t

def loadAndSplitNGon(fileName):
    dt = Filter2.loadAsChunks(fileName)
    arrays = []
    zones = I.getZones(dt)

    z = zones[0]
    cx = I.getNodeFromName2(z, 'CoordinateX')[1]
    cy = I.getNodeFromName2(z, 'CoordinateY')[1]
    cz = I.getNodeFromName2(z, 'CoordinateZ')[1]

    ngon = I.getNGonNode(z)
    ngonc = I.getNodeFromName1(ngon, 'ElementConnectivity')[1]
    ngonso = I.getNodeFromName1(ngon, 'ElementStartOffset')[1]

    nface = I.getNFaceNode(z)
    nfacec = I.getNodeFromName1(nface, 'ElementConnectivity')[1]
    nfaceso = I.getNodeFromName1(nface, 'ElementStartOffset')[1]

    fsolc = I.getNodeFromName2(z, I.__FlowSolutionCenters__)
    solc = []; solcNames = []
    if fsolc is not None:
        for f in fsolc[2]:
            if f[3] == 'DataArray_t':
                solc.append(f[1]); solcNames.append(f[0])

    fsol = I.getNodeFromName2(z, I.__FlowSolutionNodes__)
    soln = []; solNames = []
    if fsol is not None:
        for f in fsol[2]:
            if f[3] == 'DataArray_t':
                soln.append(f[1]); solNames.append(f[0])

    zonebc = I.getNodeFromType1(z, 'ZoneBC_t')
    bcs = []
    bcNames = []
    bcTypes = {}
    bcTags = {}
    if zonebc is not None:
        BCs = I.getNodesFromType1(zonebc, 'BC_t')
        for i in range(len(BCs)):
            bc = BCs[i]
            bcname = bc[0]
            bctype = I.getValue(bc)

            if bctype == 'FamilySpecified':
                fname = I.getNodeFromType(bc, 'FamilyName_t')
                fn = I.getValue(fname)
                bcTypes[bcname] = fn
            else:
                bcTypes[bcname] = bctype

            bcNames.append(bcname)
            bcTags[bcname] = i

            plist = I.getNodeFromName1(bc, 'PointList')
            bcs.append(plist[1][0])

    arrays.append([cx,cy,cz,ngonc,ngonso,nfacec,nfaceso,solc,soln,bcs])

    RES = xcore.chunk2partNGon(arrays)
    (mesh, comm_data, solc, sol, bcs, cells, faces, points) = RES

    Cmpi.barrier()

    # create zone
    zo = I.createZoneNode('{}_{}'.format(I.getName(z), Cmpi.rank), mesh)

    # add ZoneGridConnectivity
    ZGC = I.newZoneGridConnectivity(parent=zo)

    for data in comm_data:
        Name = 'Match_'+str(data[0])
        I.newGridConnectivity1to1(name=Name, donorName=str(data[0]), pointList=data[1], parent=ZGC)

    I.newUserDefinedData(name='CellLoc2Glob', value=RES[5], parent=ZGC)
    I.newUserDefinedData(name='FaceLoc2Glob', value=RES[6], parent=ZGC)
    I.newUserDefinedData(name='PointLoc2Glob', value=RES[7], parent=ZGC)

    # add solutions
    for n, name in enumerate(solNames):
        cont = I.createUniqueChild(zo, I.__FlowSolutionNodes__, 'FlowSolution_t')
        I.newDataArray(name, value=sol[n], parent=cont)

    for n, name in enumerate(solcNames):
        cont = I.createUniqueChild(zo, I.__FlowSolutionCenters__, 'FlowSolution_t')
        I._createUniqueChild(cont, 'GridLocation', 'GridLocation_t', value='CellCenter')
        I.newDataArray(name, value=solc[n], parent=cont)

    # add bcs
    for i, bc in enumerate(bcs):
        if len(bc) != 0:
            cont = I.createUniqueChild(zo, 'ZoneBC', 'ZoneBC_t')
            val = bcTypes[bcNames[i]]
            tag = bcTags[bcNames[i]]
            bcn = []
            if val not in BCType_l:
                bcn = I.newBC(name=bcNames[i], pointList=bc, family=val, parent=cont)
            else:
                bcn = I.newBC(name=bcNames[i], pointList=bc, btype=val, parent=cont)
            I.newUserDefinedData(name='Tag', value=tag, parent=bcn)

    t = C.newPyTree(['Base', zo])
    Cmpi._setProc(t, Cmpi.rank)

    # copy families
    base = I.getNodeFromName1(t, 'Base')
    families = I.getNodesFromType2(dt, 'Family_t')
    for fam in families:
        I.duptree__(fam, base)

    I._correctPyTree(t, level=7)

    return t, RES

######################################################
######################################################
# IntersectMesh functions
######################################################
######################################################
def IntersectMesh_Exit(IM):
    return xcore.IntersectMesh_Exit(IM)

def IntersectMesh_Init(t):
    zones = I.getZones(t)
    z = zones[0]
    array = C.getFields(I.__GridCoordinates__, z, api=3)[0]

    tags = I.getNodeFromName(z, 'keep')
    if tags is not None: tags = tags[1]

    return xcore.IntersectMesh_Init(array, tags)

def IntersectMesh_ExtractMesh(IM, removePeriodic=0):
    marray = xcore.IntersectMesh_ExtractMesh(IM, removePeriodic)
    zone = I.createZoneNode("Mesh", marray)
    t = C.newPyTree(["Mesh", zone])
    try: import Intersector.PyTree as XOR
    except: raise ImportError("XCore.PyTree: requires Intersector.PyTree module.")
    t = XOR.closeCells(t)
    return t

def removeIntersectingKPlanes(IM, slave_struct):
    slave_bases = I.getBases(slave_struct)

    iter = -1

    import Generator.PyTree as G
    import Transform.PyTree as T

    ts = I.newCGNSTree()

    for slave_base in slave_bases:

        iter = iter + 1

        zs = I.getZones(slave_base)

        slaves = []

        print("doing base" + str(iter))

        for z in zs:
            smesh = C.getFields(I.__GridCoordinates__, z, api=3)[0]
            slaves.append(smesh)

        new_slaves_and_tags = xcore.removeIntersectingKPlanes(IM, slaves)

        new_base = I.newCGNSBase('slave'+str(iter), 3, 3, parent=ts)

        zones = []
        for i in range(len(new_slaves_and_tags)):
            new_slave, tag = new_slaves_and_tags[i]
            zname = zs[i][0]
            zo = I.createZoneNode(zname, new_slave)
            cont = I.createUniqueChild(zo, I.__FlowSolutionNodes__, 'FlowSolution_t')
            I.newDataArray("tag", value=tag, parent=cont)
            C._convertArray2NGon(zo)
            G._close(zo)
            zones.append(zo)

        merged = T.merge(zones)
        merged = G.close(merged)
        I.addChild(new_base, merged)


    return ts

def prepareMeshesForIntersection(IM, slave):
    zs = I.getZones(slave)[0]

    s = C.getFields(I.__GridCoordinates__, zs, api=3)[0]

    tag = I.getNodeFromName2(zs, "tag")
    if tag is None:
        raise ValueError("Tag field not found in slave mesh.")
    tag = I.getValue(tag)

    keep = I.getNodeFromName(zs, 'keep')

    m, mpatch, s, spatch = xcore.prepareMeshesForIntersection(IM, s, tag)

    zmo = I.createZoneNode("M_adapted", m)
    zbcs = I.createUniqueChild(zmo, 'ZoneBC', 'ZoneBC_t')
    I.newBC(name="intersection_patch", pointList=mpatch, family='UserDefined', parent=zbcs)
    ma = C.newPyTree(["M_adapted", zmo])

    zso = I.createZoneNode("S_adapted", s)
    zbcs = I.createUniqueChild(zso, 'ZoneBC', 'ZoneBC_t')
    I.newBC(name="intersection_patch", pointList=spatch, family='UserDefined', parent=zbcs)
    sa = C.newPyTree(["S_adapted", zso])

    try: import Intersector.PyTree as XOR
    except: raise ImportError("XCore.PyTree: requires Intersector.PyTree module.")

    ma = XOR.closeCells(ma)
    sa = XOR.closeCells(sa)

    if keep is not None:
        C._cpVars(slave, 'centers:keep', sa, 'centers:keep')

    return ma, sa

def intersectMesh(master, slave):
    zm = I.getZones(master)[0]
    marr = C.getFields(I.__GridCoordinates__, zm, api=3)[0]
    mpatch = I.getNodeFromName(zm, "intersection_patch")[2][0][1]

    zs = I.getZones(slave)[0]
    keep = I.getNodeFromName(zs, 'keep')
    sarr = C.getFields(I.__GridCoordinates__, zs, api=3)[0]
    #spatch = I.getNodeFromName(zs, "intersection_patch")[2][1][1][0]
    spatch = I.getNodeFromName(zs, "intersection_patch")[2][0][1]

    marr, sarr = xcore.intersectMesh(marr, mpatch, sarr, spatch)

    zmo = I.createZoneNode("mi", marr)
    zso = I.createZoneNode("si", sarr)

    mi = C.newPyTree(["mi", zmo])
    si = C.newPyTree(["si", zso])

    try: import Intersector.PyTree as XOR
    except: raise ImportError("XCore.PyTree: requires Intersector.PyTree module.")

    mi = XOR.closeCells(mi)
    si = XOR.closeCells(si)

    if keep is not None:
        C._cpVars(slave, 'centers:keep', si, 'centers:keep')

    return mi, si

###############################################################################

def splitConnex(m):
    zones = I.getZones(m)
    if len(zones) != 1: raise ValueError('Master should be one zone.')
    zm = zones[0]
    marr = C.getFields(I.__GridCoordinates__, zm, api=3)[0]
    ptag = I.getNodeFromName(zm, 'tag')
    if ptag is None: raise ValueError('Missing point tags')
    ctag = I.getNodeFromName(zm, 'keep')
    if ctag is None: raise ValueError('Missing cell tags')
    new_arrs, new_ctags, new_ptags = xcore.split_connex(marr, ctag[1], ptag[1])
    zout = []
    for i in range(len(new_arrs)):
        z = I.createZoneNode("zone"+str(i), new_arrs[i])
        cont = I.createUniqueChild(z, I.__FlowSolutionCenters__, 'FlowSolution_t')
        I._createUniqueChild(cont, 'GridLocation', 'GridLocation_t', value='CellCenter')
        I.newDataArray('keep', value=new_ctags[i], parent=cont)
        cont = I.createUniqueChild(z, I.__FlowSolutionNodes__, 'FlowSolution_t')
        I.newDataArray('tag', value=new_ptags[i], parent=cont)
        zout.append(z)
    return zout

def icapsuleInit2():
    return xcore.icapsule_init2()

def icapsuleAdapt2(IC):
    return xcore.icapsule_adapt2(IC)

def icapsuleIntersect2(IC):
    return xcore.icapsule_intersect2(IC)

def icapsuleSetMaster(IC, m):
    zones = I.getZones(m)
    if len(zones) != 1: raise ValueError('Master should be one zone.')
    zm = zones[0]
    marr = C.getFields(I.__GridCoordinates__, zm, api=3)[0]
    ctag = I.getNodeFromName(zm, 'keep')
    if ctag is None: raise ValueError('Missing cell tags')
    return xcore.icapsule_set_master(IC, marr, ctag[1])

def icapsuleSetSlaves(IC, slaves):
    sarrs = []
    ptags = []
    ctags = []

    for slave in slaves:

        bases = I.getBases(slave)

        for base in bases:
            zones = I.getZones(base)
            for zone in zones:
                sarr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
                sarrs.append(sarr)
                ptag = I.getNodeFromName(zone, 'tag')
                if ptag is None: raise ValueError('Missing point tags.')
                ptags.append(ptag[1])
                ctag = I.getNodeFromName(zone, 'keep')
                if ctag is None: raise ValueError('Missing cell tags.')
                ctags.append(ctag[1])

    return xcore.icapsule_set_slaves(IC, sarrs, ptags, ctags)

# Extract the master mesh corresponding to the capsule
def icapsuleExtractMaster(IC):
    marr, ctag = xcore.icapsule_extract_master(IC)
    zm = I.createZoneNode("master", marr)

    cont = I.createUniqueChild(zm, I.__FlowSolutionCenters__, 'FlowSolution_t')
    I._createUniqueChild(cont, 'GridLocation', 'GridLocation_t', value='CellCenter')
    I.newDataArray("keep", value=ctag, parent=cont)
    return zm

# Extract the slave mesh corresponding to the capsule
def icapsuleExtractSlaves(IC):
    sarrs, ctags = xcore.icapsule_extract_slaves(IC)
    assert(len(sarrs) == len(ctags))
    zones = []
    for i in range(len(sarrs)):
        zs = I.createZoneNode("slave_"+str(i), sarrs[i])

        cont = I.createUniqueChild(zs, I.__FlowSolutionCenters__, 'FlowSolution_t')
        I._createUniqueChild(cont, 'GridLocation', 'GridLocation_t', value='CellCenter')
        I.newDataArray("keep", value=ctags[i], parent=cont)

        zones.append(zs)
    return zones

# Triangulates the external quads of a mesh
# and modifies the corresponding BC PointList
def triangulateSkin(m):
    m_copy = I.copyRef(m)
    _triangulateSkin(m_copy)
    return m_copy

def _triangulateSkin(m):
    zones = I.getZones(m)
    for zone in zones:
        marr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
        zbc = I.getNodeFromType1(zone, 'ZoneBC_t')
        ptlists = []
        if zbc is not None:
            bcs = I.getNodesFromType1(zbc, 'BC_t')
            for bc in bcs:
                ptlists.append(I.getNodeFromName1(bc, 'PointList')[1][0])
        m_out, ptlists_out = xcore.triangulate_skin(marr, ptlists)
        C.setFields([m_out], zone, 'nodes')
        if zbc is not None:
            bcs = I.getNodesFromType1(zbc, 'BC_t')
            for j, bc in enumerate(bcs):
                ptlist = I.getNodeFromName1(bc, 'PointList')
                ptlist[1] = ptlists_out[j]

    return None
