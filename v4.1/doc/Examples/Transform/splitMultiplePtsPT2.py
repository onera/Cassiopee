# split des blocs en raccord multiple 2D
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

z0 = G.cart((0.,0.,0.),(0.1,0.1,1.),(10,10,1));
z1 = T.subzone(z0,(1,1,1),(5,10,1)); z1[0] = 'cart1'
z1 = C.addBC2Zone(z1, 'wall1','BCWall','jmin')
z2 = T.subzone(z0,(5,1,1),(10,5,1)); z2[0] = 'cart2'
z3 = T.subzone(z0,(5,5,1),(10,10,1)); z3[0] = 'cart3'
t = C.newPyTree(['Base',2])
t[2][1][2] += [z1, z2, z3]
t = T.splitMultiplePts(t)
C.convertPyTree2File(t, "out.cgns", "bin_cgns")
