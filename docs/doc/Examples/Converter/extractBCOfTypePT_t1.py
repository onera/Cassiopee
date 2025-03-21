# - extractBCOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

# Sur une zone structure + champ en noeuds + champ en centres
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
Z = C.extractBCOfType(a, 'BCWall')
t = C.newPyTree(['Base',3,'Skin',2])
t[2][1][2].append(a); t[2][2][2] += Z
test.testT(t,1)

# Sur un arbre de zones structurees
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
a = C.addBC2Zone(a,'wall','BCWall','imin')
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
b = G.cart((0,0,0), (0.1,0.1,0.1), (10,30,10))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
t = C.newPyTree(['Base',3]); t[2][1][2] += [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
Z = C.extractBCOfType(t, 'BCWall')
t = C.newPyTree(['Skin',2]); t[2][1][2] += Z
test.testT(t, 2)

# Sur une liste de zones structurees
a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
b = G.cart((0,0,0), (0.1,0.1,0.1), (10,30,10))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'jmin')
Z = C.extractBCOfType([a,b], 'BCWall')
t = C.newPyTree(['Skin',2]); t[2][1][2] += Z
test.testT(t, 3)

# BC referencees par des familles
a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
a = C.addBC2Zone(a, 'walla', 'FamilySpecified:CARTER', 'imin')
a = C.addBC2Zone(a, 'nref', 'FamilySpecified:LOIN', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a]
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER', bndType='BCWall')
t[2][1] = C.addFamily2Base(t[2][1], 'LOIN', bndType='BCFarfield')
zones = C.extractBCOfType(t, 'BCWall')
t = C.newPyTree(['Base']); t[2][1][2] += zones
test.testT(t, 4)

# BC sur un NGON
a = G.cartNGon( (0,0,0), (1,1,1), (10,10,10) )
subzone = G.cartNGon( (0,0,0), (1,1,1), (10,10,1) )
hook = C.createHook(a, function='faceCenters')
ids = C.identifyElements(hook, subzone)
C._addBC2Zone(a, 'wall', 'BCWall', faceList=ids)
ext = C.extractBCOfType(a, 'BCWall')
test.testT(ext, 5)

# BC sur un basic elements
a = G.cartHexa( (0,0,0), (1,1,1), (5,5,5) )
subzone = G.cartHexa( (0,0,0), (1,1,1), (5,5,1) )
hook = C.createHook(a, function='faceCenters')
ids = C.identifyElements(hook, subzone)
a = C.addBC2Zone(a, 'wall', 'BCWall', faceList=ids)
ext = C.extractBCOfType(a, 'BCWall')
test.testT(ext, 6)

# BC sur un BE + connectivite QUAD de boundary
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartHexa((0,0,0), (1,1,1), (10,10,1))
C._addBC2Zone(a, 'wall', 'BCWall', subzone=b)
bc = Internal.getNodeFromType(a,'BC_t')
ER = Internal.getNodeFromName(bc,'ElementRange')
ER = Internal.getValue(ER)[0]
ermin = ER[0]
ermax = ER[1]
Internal._rmNodesFromName(bc,'ElementRange')
import numpy
nfaces = ermax-ermin+1
r = numpy.zeros(nfaces, dtype=Internal.E_NpyInt)
for i in range(nfaces): r[i]=ermin+i
Internal._createChild(bc, 'GridLocation', 'GridLocation_t', value='FaceCenter')
r = r.reshape((1,r.size), order='F')
bc[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
ext = C.extractBCOfType(a, 'BCWall')
test.testT(ext, 7)
