import Converter.Filter2 as Filter2
import Converter.Mpi as Cmpi
import Converter.Internal as I
import Converter.PyTree as C
import XCore.xcore

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