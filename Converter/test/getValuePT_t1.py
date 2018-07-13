# - getValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Array structure
Ni = 40; Nj = 50; Nk = 20
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
# Les variables contenues dans a (x,y,z) au point (10,1,1)
val1 = C.getValue(a, 'CoordinateX', (10,1,1))
val2 = C.getValue(a, 'CoordinateX', 9) # C'est le meme point
val3 = C.getValue(a, 'nodes:CoordinateX', 9) # C'est le meme point
val4 = C.getValue(a, 'GridCoordinates', 9) # retourne (x,y,z)
test.testO([val1, val2, val3, val4], 1)

# Array structure + champ en centres
a = G.cart((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
a = C.initVars(a, 'centers:Density', 1)
val = C.getValue(a, 'centers:Density', (1,1,1) )
test.testO([val], 11)

# Array non structure
Ni = 40; Nj = 50; Nk = 20
a = G.cartTetra((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
val = C.getValue( a, 'CoordinateX', 9 )
test.testO(val, 2)

# Array non structure + champs en centre
a = G.cartTetra((0,0,0), (1./(Ni-1), 0.5/(Nj-1),1./(Nk-1)), (Ni,Nj,Nk))
a = C.initVars(a, 'centers:Density', 1)
val = C.getValue(a, 'centers:Density', 10 )
test.testO([val], 21)

# test sur un arbre
b = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 2))
b = C.addBC2Zone(b, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t = C.initVars(t,'F',1.); t = C.initVars(t,'centers:G',2.)
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)
val = C.getValue( a, 'CoordinateX', 9 )
test.testO(val, 3)
