import Converter.Filter2 as Filter2
import Converter.Mpi as Cmpi
import Converter.Internal as I
import Converter.PyTree as C
import XCore.xcore

def chunk2part(dt):
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

  out = []

  for cn in cns:
    name, stride = I.eltNo2EltName(cn[1][0])
    arr = I.getNodeFromName1(cn, 'ElementConnectivity')[1]
    mesh = XCore.xcore.chunk2part_elt(XYZ, name, stride, arr)
    out.append(mesh)

  return out
  