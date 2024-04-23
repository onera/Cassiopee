import Converter.Filter2 as Filter2
import Converter.Mpi as Cmpi
import Converter.Internal as I
import Converter.PyTree as C
from . import xcore

# Returns for each zone, exchanged fields
def exchangeFields(t, fldNames):
    zones = I.getZones(t)
    rfields = []
    for zone in zones:
        arr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
        pe = I.getNodeFromName(zone, 'ParentElements')
        if pe == None: raise ValueError('ParentElements not found.')
        fsolc = I.getNodeFromName2(zone, I.__FlowSolutionCenters__)
        if fsolc == None: raise ValueError('FlowSolutionCenters not found.')
        flds = []
        for fldName in fldNames:
            fld = I.getNodeFromName2(fsolc, fldName)
            if fld == None: raise ValueError(fldName, 'not found.')
            flds.append(fld[1])
        zgc = I.getNodeFromType(zone, 'ZoneGridConnectivity_t')
        if zgc == None: raise ValueError('ZoneGridConnectivity not found')
        comms = I.getNodesFromType(zgc, 'GridConnectivity1to1_t')
        if comms == None: raise ValueError('GridConnectivity1to1 not found')
        comm_list = []
        for comm in comms:
            nei_proc = int(I.getValue(comm))
            ptlist = I.getNodeFromName(comm, 'PointList')[1]
            comm_list.append([nei_proc, ptlist])
        rfields.append(xcore.exchangeFields(arr, pe[1], flds, comm_list))
    return rfields

def initAdaptTree(t):
  zones = I.getZones(t)
  adaptTrees = []
  for z in zones:
    fc = C.getFields(I.__GridCoordinates__, z, api=3)[0]
    if fc != []:
      adaptTrees.append(xcore.initAdaptTree(fc))
    else:
      adaptTrees.append(None)
  return adaptTrees

def loadAndSplitElt(fileName):
  dt = Filter2.loadAsChunks(fileName)
  zones = I.getZones(dt)
  
  if len(zones) > 1:
    raise TypeError("chunk2part: one zone only.")

  z = zones[0]

  XYZ = []

  cx = I.getNodeFromName2(z, 'CoordinateX')[1]
  cy = I.getNodeFromName2(z, 'CoordinateY')[1]
  cz = I.getNodeFromName2(z, 'CoordinateZ')[1]

  XYZ.append(cx)
  XYZ.append(cy)
  XYZ.append(cz)

  cns = I.getNodesFromType(z, 'Elements_t')

  chunks = []

  for cn in cns:
    name, stride = I.eltNo2EltName(cn[1][0])
    arr = I.getNodeFromName1(cn, 'ElementConnectivity')[1]
    chunks.append([name, stride, arr])

  parts = xcore.chunk2partElt(XYZ, chunks)

  zones = []

  for i in range(len(parts)):
    z = I.createZoneNode('Zone' + '%d'%Cmpi.rank + '_%d'%i, parts[i])
    zones.append(z)

  t = C.newPyTree(['Base', zones])

  Cmpi._setProc(t, Cmpi.rank)

  return t

BCType_l = set(I.KNOWNBCS)

def loadAndSplitNGon(fileName):
  dt = Filter2.loadAsChunks(fileName)
  arrays = []
  zones = I.getZones(dt)

  z = zones[0]
  cx = I.getNodeFromName2(z, 'CoordinateX')[1]
  cy = I.getNodeFromName2(z, 'CoordinateY')[1]
  cz = I.getNodeFromName2(z, 'CoordinateZ')[1]

  ngon = I.getNodeFromName2(z, 'NGonElements')
  ngonc = I.getNodeFromName1(ngon, 'ElementConnectivity')[1]
  ngonso = I.getNodeFromName1(ngon, 'ElementStartOffset')[1]

  nface = I.getNodeFromName2(z, 'NFaceElements')
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

  zonebc = I.getNodeFromType(z, 'ZoneBC_t')
  bcs = []
  bcNames = []
  bcTypes = {}
  if zonebc is not None:
    BCs = I.getNodesFromType(zonebc, 'BC_t')
    for bc in BCs:
      bcname = bc[0]
      bctype = I.getValue(bc)

      if bctype == 'FamilySpecified':
        fname = I.getNodeFromType(bc, 'FamilyName_t')
        fn = I.getValue(fname)
        bcTypes[bcname] = fn
      else:
        bcTypes[bcname] = bctype

      bcNames.append(bcname)

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
    I._createUniqueChild(cont, 'GridLocation', 'GridLocation_t', value='CellCenter', )
    I.newDataArray(name, value=solc[n], parent=cont)
  
  for i in range(len(bcs)):
    if len(bcs[i]) != 0:
      cont = I.createUniqueChild(zo, 'ZoneBC', 'ZoneBC_t')
      val = bcTypes[bcNames[i]]
      if val not in BCType_l:
        I.newBC(name=bcNames[i], pointList=bcs[i], family=val, parent=cont)
      else:
        I.newBC(name=bcNames[i], pointList=bcs[i], btype=val, parent=cont)

  t = C.newPyTree(['Base', zo])
  Cmpi._setProc(t, Cmpi.rank)
  
  # copy families
  base = I.getNodeFromName1(t, 'Base')
  families = I.getNodesFromType2(dt, 'Family_t')
  for fam in families:
    I.duptree__(fam, base)
  
  I._correctPyTree(t, level=7)

  return t, RES

def _adaptMeshDir(h, l, fld):
    zone = I.getZones(l)[0]
    arr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
    xcore.adaptMeshDir(h, arr, fld)
    return None

def _adaptMeshSeq(h, fld, fv=None):
    if isinstance(h, list):
        if len(h) != len(fld): raise ValueError('mesh hooks and fields not the same size')
        for i in range(len(h)):
            xcore.adaptMeshSeq(h[i], fld[i], fv)
    else:
        xcore.adaptMeshSeq(h, fld, fv)
    return None

def extractLeafMesh(h):
    if isinstance(h, list):
        leaves = []
        for i in range(len(h)):
            m = xcore.extractLeafMesh(h[i])
            leaves.append(I.createZoneNode('Leaves' + '%d'%i, m))
        T = C.newPyTree(['Base', leaves])
        return T
    else:
        m = xcore.extractLeafMesh(h)
        leaf = I.createZoneNode('Leaves', m)
        T = C.newPyTree(['Base', leaf])
        return T

def createAdaptMesh(t):
    I._adaptNGon32NGon4(t)
    zones = I.getZones(t)
    AMs = []
    for z in zones:
        fc = C.getFields(I.__GridCoordinates__, z, api=3)[0]
        if fc != []:
            AMs.append(xcore.createAdaptMesh(fc))
    return AMs

############################################################################

def CreateAdaptMesh(t, own, nei, comm, Tr, Tu, eps, hmin, hmax, unrefine=False,
  mode_2D=None, gcells=None, gfaces=None, gpoints=None):
  zones = I.getZones(t)
  z = zones[0]
  fc = C.getFields(I.__GridCoordinates__, z, api=3)[0]

  '''
  comm_list = []
  zgc = I.getNodeFromType(z, 'ZoneGridConnectivity_t')
  if zgc is not None:
    comm_data = I.getNodesFromType(zgc, 'GridConnectivity1to1_t')
    
    for data in comm_data:
      nei_proc = int(I.getValue(data))
      pfaces = I.getNodeFromName(data, 'PointList')[1]
      comm_list.append([nei_proc, pfaces])
  '''
  
  bcs = []
  zonebc = I.getNodeFromType(z, 'ZoneBC_t')
  if zonebc is not None:
    zbc = I.getNodesFromType(zonebc, 'BC_t')
  
    for bc in zbc:
      plist = I.getNodeFromName(bc, 'PointList')
      #name = I.getNodeFromName(bc, 'FamilyName')
      name = bc[0]
      #name = str(I.getValue(name))
      bcs.append([plist[1], name])
    
  return xcore.CreateAdaptMesh(fc, own, nei, comm, bcs, Tr, Tu, eps,
    hmin, hmax, unrefine, mode_2D, gcells, gfaces, gpoints)

def AdaptMesh(AM):
  return xcore.AdaptMesh(AM)

def computeGradient(AM, field, cx, cy, cz, own, nei):
  return xcore.computeGradient(AM, field, cx, cy, cz, own, nei)

def computeHessian(AM, field, grad, cx, cy, cz, own, nei):
  return xcore.computeHessian(AM, field, grad, cx, cy, cz, own, nei)

def hessianToMetric(H, hmin, hmax, eps):
  return xcore.hessianToMetric(H, hmin, hmax, eps)

def _makeRefDataFromGradAndHess(AM, f, g, h):
  return xcore._makeRefDataFromGradAndHess(AM, f, g, h)

def _prepareMeshForAdaptation(t):
    zones = I.getZones(t)
    zone = zones[0]
    arr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
    return xcore._prepareMeshForAdaptation(arr)

def ExtractLeafMesh(t, conformize=1):
  return xcore.ExtractLeafMesh(t, conformize)

def _assignRefDataToAM(AM, REF):
  return xcore._assignRefDataToAM(AM, REF)

def extractBoundaryMesh(AM, mode):
  return xcore.extractBoundaryMesh(AM, mode)


######################################################

def intersectMesh(master, slave, patch_name):
  return removeIntersectingKPlanes(master, slave, patch_name)

def removeIntersectingKPlanes(master, slave, patch_name):
  zm = I.getZones(master)[0]
  zs = I.getZones(slave)[0]

  master = C.getFields(I.__GridCoordinates__, zm, api=3)[0]
  slave = C.getFields(I.__GridCoordinates__, zs, api=3)[0]

  patch = I.getNodeFromName(zm, patch_name)
  if patch is None:
    raise ValueError(patch_name + "not found.")
  
  faces = I.getNodeFromName(patch, "PointList")
  faces = I.getValue(faces)[0]

  mesh, tag = xcore.removeIntersectingKPlanes(master, slave, faces)

  zo = I.createZoneNode("struct", mesh)
  cont = I.createUniqueChild(zo, I.__FlowSolutionNodes__, 'FlowSolution_t')
  I.newDataArray("tag", value=tag, parent=cont)

  t = C.newPyTree(["Base", zo])
  return t

def intersectSurf(master, slave, patch_name):
  zm = I.getZones(master)[0]
  zs = I.getZones(slave)[0]

  m = C.getFields(I.__GridCoordinates__, zm, api=3)[0]
  s = C.getFields(I.__GridCoordinates__, zs, api=3)[0]

  patch = I.getNodeFromName(zm, patch_name)
  if patch is None:
    raise ValueError(patch_name + "not found.")
  
  faces = I.getNodeFromName(patch, "PointList")
  faces = I.getValue(faces)[0]

  tag = I.getNodeFromName2(zs, "tag")
  if tag is None:
    raise ValueError("Tag field not found in slave mesh.")
  tag = I.getValue(tag)

  mesh = xcore.intersectSurf(m, s, faces, tag)

  return mesh
  

def sweep(m, s):
    zm = I.getZones(m)[0]
    zs = I.getZones(s)[0]

    master = C.getFields(I.__GridCoordinates__, zm, api=3)[0]
    slave = C.getFields(I.__GridCoordinates__, zs, api=3)[0]

    return xcore.sweep(master, slave)
