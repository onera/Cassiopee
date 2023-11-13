import Converter.Filter2 as Filter2
import Converter.Mpi as Cmpi
import Converter.Internal as I
import Converter.PyTree as C
import XCore.xcore

def exchangeFields(t, fldnames):
    zones = I.getZones(t)
    for zone in zones:
        arr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
        pe = I.getNodeFromName(zone, 'ParentElements')
        if pe == None: raise ValueError('ParentElements not found.')
        fsolc = I.getNodeFromName2(zone, I.__FlowSolutionCenters__)
        if fsolc == None: raise ValueError('FlowSolutionCenters not found.')
        flds = []
        for fldname in fldnames:
            fld = I.getNodeFromName2(fsolc, fldname)
            if fld == None: raise ValueError(fldname + 'not found.')
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
        rfields = XCore.xcore.exchangeFields(arr, pe[1], flds, comm_list)
    
    return t

def initAdaptTree(t):
  zones = I.getZones(t)
  adaptTrees = []
  for z in zones:
    fc = C.getFields(I.__GridCoordinates__, z, api=3)[0]
    if fc != []:
      adaptTrees.append(XCore.xcore.initAdaptTree(fc))
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

  parts = XCore.xcore.chunk2partElt(XYZ, chunks)

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
  familyNames = {}
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

  RES = XCore.xcore.chunk2partNGon(arrays)

  mesh = RES[0]
  comm_data = RES[1]
  solc = RES[2]
  sol = RES[3]
  bcs = RES[4]
  cells = RES[5]
  faces = RES[6]
  points = RES[7]

  Cmpi.barrier()

  # create zone
  zo = I.createZoneNode('Zone_' + '%d'%Cmpi.rank, mesh)

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
  I._correctPyTree(t, level=7)

  return t, RES

def _adaptMeshDir(h, l, fld):
    zone = I.getZones(l)[0]
    arr = C.getFields(I.__GridCoordinates__, zone, api=3)[0]
    XCore.xcore.adaptMeshDir(h, arr, fld)
    return None

def _adaptMeshSeq(h, fld, fv=None):
    if isinstance(h, list):
        if len(h) != len(fld): raise ValueError('mesh hooks and fields not the same size')
        for i in range(len(h)):
            XCore.xcore.adaptMeshSeq(h[i], fld[i], fv)
    else:
        XCore.xcore.adaptMeshSeq(h, fld, fv)
    return None

def extractLeafMesh(h):
    if isinstance(h, list):
        leaves = []
        for i in range(len(h)):
            m = XCore.xcore.extractLeafMesh(h[i])
            leaves.append(I.createZoneNode('Leaves' + '%d'%i, m))
        T = C.newPyTree(['Base', leaves])
        return T
    else:
        m = XCore.xcore.extractLeafMesh(h)
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
            AMs.append(XCore.xcore.createAdaptMesh(fc))
    return AMs
