# - hasBCDataSet (pyTree) -
import Converter.Internal as Internal
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# -- BCDataSet
# On a zone
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.addBC2Zone(a, 'wall', 'BCWall', 'imin')
res = Internal.hasBCDataSets(a)
t = C.newPyTree(['Base', a])
test.testO(res, 1)

b = Internal.getNodeFromName2(a, 'wall')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='FaceCenter', parent=b)
d = Internal.newBCData('BCNeumann', parent=d)
d = Internal.newDataArray('Density', value=10*[1.], parent=d)
res = Internal.hasBCDataSets(a)
test.testO(res, 2)

res = Internal.hasBCDataSets(a, loc='Vertex')
test.testO(res, 3)

# On a PyTree
res = Internal.hasBCDataSets(t, loc='FaceCenter')
test.testO(res, 4)

res = Internal.hasBCDataSets(t, loc='Vertex')
test.testO(res, 5)

# On a list of zones
a = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((5,0.9,0), (0.01,0.01,1.), (20,20,2)); b[0] = 'cart2'
zones = [a, b]
C._addBC2Zone(zones, 'wall1', 'BCWallInviscid', 'imin')
res = Internal.hasBCDataSets(zones, loc='FaceCenter')
test.testO(res, 6)

b = Internal.getNodeFromName(zones, 'wall1')
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                          gridLocation='Vertex', parent=b)
d = Internal.newBCData(parent=d)
d = Internal.newDataArray('Density', value=20*[1.], parent=d)

res = Internal.hasBCDataSets(zones, loc='FaceCenter')
test.testO(res, 7)

res = Internal.hasBCDataSets(zones, loc='Vertex')
test.testO(res, 8)

# -- SubRegion
"""
n = Internal.newZoneSubRegion(name='SubRegionBC', bcName='BC', gridLocation='FaceCenter'); Internal.printTree(n)
print("---- SubRegion")
a = G.cart((0,0,0), (1,1,1), (10,10,1))
n = Internal.newZoneSubRegion(name='SubRegionBC', bcName='BC', gridLocation='FaceCenter')
res = Internal.hasBCDataSets(a)
print(n)
Internal.printTree(a)

res = Internal.hasBCDataSets(a, loc='FaceCenter')
print(res)
res = Internal.hasBCDataSets(a, loc='Vertex')
print(res)
"""
