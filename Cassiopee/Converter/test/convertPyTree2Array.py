# - convertPyTree2Array (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
a = C.initVars(a,'F',1.); a = C.initVars(a,'centers:G',2.)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'imin')
b = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(11,11,11)); b[0] = 'cartHexa'
t = C.newPyTree(['Base']); t[2][1][2] = t[2][1][2] + [a,b]
t[2][1] = C.addState(t[2][1], 'Mach', 0.6)

zones = C.convertPyTree2ZoneNames(t)

# 1 - Get one field
arrays = []
for i in zones:
    a = C.convertPyTree2Array(i+"/GridCoordinates/CoordinateX", t)
    b = C.convertPyTree2Array(i+"/GridCoordinates/CoordinateY", t)
    c = C.convertPyTree2Array(i+"/GridCoordinates/CoordinateZ", t)
    x = Converter.addVars([a,b,c])
    arrays.append(x)
Converter.convertArrays2File(arrays, "out.plt")

# 2 - Get a global field
arrays = []
for i in zones:
    a = C.convertPyTree2Array(i+"/GridCoordinates", t)
    arrays.append(a)
Converter.convertArrays2File(arrays, "out2.plt")

# 3 - Get the flow solution
a = C.convertPyTree2Array(zones[0]+'/FlowSolution', t)
print(a)
