import KCore

import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy

ns = 200
a = G.cart((0,0,0), (1,1,1), (ns,ns,ns))
# Get an array2 from a
x1 = Internal.getNodeFromName2(a, 'CoordinateX')
x2 = Internal.getNodeFromName2(a, 'CoordinateY')
x3 = Internal.getNodeFromName2(a, 'CoordinateZ')
b = ['x,y,z', [x1[1],x2[1],x3[1]], ns, ns, ns]
#KCore.tester(b)

a = G.cartHexa((0,0,0), (1,1,1), (ns,ns,ns))
# Get an array2 from a
x1 = Internal.getNodeFromName2(a, 'CoordinateX')
x2 = Internal.getNodeFromName2(a, 'CoordinateY')
x3 = Internal.getNodeFromName2(a, 'CoordinateZ')
c1 = Internal.getNodeFromName2(a, 'ElementConnectivity')
c1 = c1[1].reshape((8,(ns-1)*(ns-1)*(ns-1))) # cool
b = ['x,y,z', [x1[1],x2[1],x3[1]], [c1], 'HEXA']
KCore.tester(b)

a = G.cartNGon((0,0,0), (1,1,1), (ns,ns,ns))
C.convertPyTree2File(a, 'out.cgns')
# Get an array2 from a
x1 = Internal.getNodeFromName2(a, 'CoordinateX')
x2 = Internal.getNodeFromName2(a, 'CoordinateY')
x3 = Internal.getNodeFromName2(a, 'CoordinateZ')
c1 = Internal.getNodesFromName2(a, 'ElementConnectivity')
b = ['x,y,z', [x1[1],x2[1],x3[1]], [c1], 'HEXA']
KCore.tester(b)


