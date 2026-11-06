# - enforcePlusZ (pyTree) -
import Generator.PyTree as G
import KCore.test as test
import Converter.PyTree as C

# Distribution
Ni = 50; Nj = 2; Nk = 50
a = G.cart((0,0,0), (1./(Ni-1), 1., 0.5/(Nk-1)), (Ni,Nj,Nk))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2])
a = C.addBC2Zone(a, 'match2', 'BCMatch', 'imax', a, 'imin', [1,2])
a = C.addBC2Zone(a, 'wall2', 'BCWall', 'jmax')
a = C.addVars(a, 'Density'); a = C.addVars(a, 'centers:cellN')
b = G.enforcePlusZ(a, 1.e-3, 10, 20)
test.testT(b, 1)
