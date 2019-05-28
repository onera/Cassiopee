# -- internal module tests --

import Converter.Internal as I
import Converter.PyTree as C

t = C.convertFile2PyTree('SquaredNozzle-06-R.cgns', 'bin_cgns')

# zone name server
name = C.getZoneName('cart'); print(name)
name = C.getZoneName('cart'); print(name)
name = C.getZoneName('cart'); print(name)

# isTopTree
print(I.isTopTree(t), 'is True.')
print(I.isTopTree(t[2][1][2][6]), 'is False.')

# isStdNode
print(I.isStdNode(t), 'is 1.')
print(I.isStdNode(t[2][1][2][6]), 'is -1.')
print(I.isStdNode(t[2][1][2]), 'is 0.')

# getNodeFromPath
a = I.getNodeFromPath(t, 'SquaredNozzle/Zone-001/GridCoordinates'); #print(a)
a = I.getNodeFromPath(a, 'CoordinateX'); #print(a)

# getNodesFromType
res = I.getNodesFromType(t, 'GasModel_t'); #print(res)

# getNodesFromName
res = I.getNodesFromName(t[2][1][2][6], 'ZoneBC'); #print(res)

# getParentOfNode
p = I.getParentOfNode(t[2][1][2][6], res[0]); #print(p)

# eltNo2EltName et eltName2EltNo
print(I.eltNo2EltName(20)); print(I.eltName2EltNo('MIXED'))

# Retourne le nom CGNS des variables
print(I.getCGNSName('x'))

# createZoneNode from array nodes + array centers
import Generator
a = Generator.cart((0,0,0), (1,1,1), (10,10,10))
z = I.createZoneNode('Cart', a); t[2][1][2].append(z)
a = Generator.cartTetra((0,0,0), (1,1,1), (10,10,10))
z = I.createZoneNode('Tetra', a); t[2][1][2].append(z)

# getZoneNodeDim
z = t[2][1][2][18]
print(I.getZoneDim(z))
z = t[2][1][2][19]
print(I.getZoneDim(z))

# Get array from dataNode
z = t[2][1][2][18]
node = I.getNodesFromName(z, 'CoordinateX')[0]
connects = []
dim = I.getZoneDim(z)
a = I.convertDataNode2Array(node, dim, connects); #print(a)

# Get arrays from container name in zone, base or tree
a = C.getFields('GridCoordinates', z); #print(a)

# Get all standard fields in a zone, base or tree
a = C.getAllFields(z, 'nodes'); #print(a)

# Copy ref
t2 = I.copyRef(t); #print(t2)
