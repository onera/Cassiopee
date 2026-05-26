# - addBC2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# - Non structure a elements basiques -

# Ajoute la connectivite BC + la BC par element range
a = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,10))
b = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,1))
a = C.mergeConnectivity(a, b, boundary=1)
C._addBC2Zone(a, 'wall', 'BCWall', elementRange=[730,810])
test.testT(a, 1)

# Ajout par une subzone
a = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,10))
b = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,1))
a = C.addBC2Zone(a, 'wall', 'BCWall', subzone=b)
test.testT(a, 2)
# Match : element range
a = G.cartHexa((0,0,0), (0.1,0.1,0.1), (11,11,11))
b = G.cartHexa((1,0,0), (0.1,0.1,0.1), (1,11,11))
a = C.mergeConnectivity(a, b, boundary=1)
a2 = G.cartHexa((1,0,0), (0.1,0.1,0.1), (11,11,11))
a2 = C.mergeConnectivity(a2, b, boundary=1)
a = C.addBC2Zone(a, 'match', 'BCMatch', elementRange=[1001,1100],zoneDonor=a2, elementRangeDonor=[1001,1100])
a2 = C.addBC2Zone(a2, 'match', 'BCMatch', elementRange=[1001,1100],zoneDonor=a, elementRangeDonor=[1001,1100])
t = C.newPyTree(['Base',a,a2])
test.testT(t, 3)

# Match : par subzone
a1 = G.cartHexa((0,0,0), (0.1,0.1,0.1), (11,11,11))
a2 = G.cartHexa((1,0,0), (0.1,0.1,0.1), (11,11,11))
sz = G.cartHexa((1,0,0), (0.1,0.1,0.1), (1,11,11))
a1 = C.mergeConnectivity(a1, sz, boundary=1)
a2 = C.mergeConnectivity(a2, sz, boundary=1)
a1 = C.addBC2Zone(a1, 'match', 'BCMatch', subzone=sz,zoneDonor=a2)
a2 = C.addBC2Zone(a2, 'match', 'BCMatch', subzone=sz,zoneDonor=a1)
t = C.newPyTree(['Base',a1,a2])
test.testT(t, 4)
# Match : par subzone sans mergeConnectivity prealable
a1 = G.cartHexa((0,0,0), (0.1,0.1,0.1), (11,11,11))
a2 = G.cartHexa((1,0,0), (0.1,0.1,0.1), (11,11,11))
sz = G.cartHexa((1,0,0), (0.1,0.1,0.1), (1,11,11))
a1 = C.addBC2Zone(a1, 'match', 'BCMatch', subzone=sz,zoneDonor=a2)
a2 = C.addBC2Zone(a2, 'match', 'BCMatch', subzone=sz,zoneDonor=a1)
t = C.newPyTree(['Base',a1,a2])
test.testT(t, 5)


# ajout par une subzone
a = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,10))
b = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,1))
a = C.addBC2Zone(a, 'stage', 'BCStage', subzone=b)
test.testT(a,6)

# BCStage by a family name
a = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,10))
b = G.cartHexa((2,0,0), (0.1,0.1,1), (10,10,1))
a = C.addBC2Zone(a, 'stage', 'FamilySpecified:StageMxpl', subzone=b)
test.testT(a,7)
