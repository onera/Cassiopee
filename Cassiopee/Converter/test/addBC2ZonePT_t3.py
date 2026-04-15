# - addBC2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# - NGons -

# Ajout par une faceList
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=[1,2])
test.testT(a, 1)

# ajout par une subzone
b = G.cartNGon((2,0,0), (0.1,0.1,1), (10,1,2))
a = C.addBC2Zone(a, 'wall', 'BCWall', subzone=b)
test.testT(a, 2)

# raccords match - ajout par une subzone
a = G.cartNGon((0,0,0), (0.1,0.1,1), (11,11,2))
faceMatch = G.cartNGon((1,0,0), (0.1,0.1,1), (1,11,2))
b = G.cartNGon((1,0,0), (0.1,0.1,1), (11,11,2))
a = C.addBC2Zone(a, 'match', 'BCMatch', subzone=faceMatch, zoneDonor=b)
b = C.addBC2Zone(b, 'match', 'BCMatch', subzone=faceMatch, zoneDonor=a)
t = C.newPyTree(['Base',a,b,faceMatch])
test.testT(t, 3)

# raccords match - ajout par une faceList
a = G.cartNGon((0,0,0), (0.1,0.1,1), (11,11,2))
b = G.cartNGon((1,0,0), (0.1,0.1,1), (11,11,2))
faceLista = [11,22,33,44,55,66,77,88,99,110]
faceListb = [1,12,23,34,45,56,67,78,89,100]
a = C.addBC2Zone(a, 'match', 'BCMatch', faceList=faceLista, zoneDonor=b, faceListDonor=faceListb)
b = C.addBC2Zone(b, 'match', 'BCMatch', faceList=faceListb, zoneDonor=a, faceListDonor=faceLista)
t = C.newPyTree(['Base',a,b])
test.testT(t,4)

# ajout par une subzone
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,10))
b = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,1))
a = C.addBC2Zone(a, 'stage', 'BCStage', subzone=b)
test.testT(a,5)

# BCStage by a family name
a = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,10))
b = G.cartNGon((2,0,0), (0.1,0.1,1), (10,10,1))
a = C.addBC2Zone(a, 'stage', 'FamilySpecified:StageMxpl', subzone=b)
test.testT(a,6)
