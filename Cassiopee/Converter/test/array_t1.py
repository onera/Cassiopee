# - array (array) -
import Converter as C
import KCore.test as test

# Structure - array1
a = C.array('x,y,z,cellN,cellNF', 12, 9, 12) # 3d
b = C.array('x,y,z,cellN,cellNF', 12, 9, 1)  # 2d
c = C.array('x,y,z,cellN,cellNF', 12, 1, 1)  # 1d
test.testA([a,b,c], 1)

# Structure - array2
a = C.array('x,y,z,cellN,cellNF', 12, 9, 12, api=2) # 3d
b = C.array('x,y,z,cellN,cellNF', 12, 9, 1, api=2)  # 2d
c = C.array('x,y,z,cellN,cellNF', 12, 1, 1, api=2)  # 1d
#test.testA([a,b,c], 3)

# Non structure - array1
a = C.array('x,y,z,cellN,cellNF', 12, 9, 'HEXA')
b = C.array('x,y,z,cellN,cellNF', 12, 9, 'TETRA')
c = C.array('x,y,z,cellN,cellNF', 12, 9, 'PENTA')
d = C.array('x,y,z,cellN,cellNF', 12, 9, 'PYRA')
e = C.array('x,y,z,cellN,cellNF', 12, 9, 'QUAD')
f = C.array('x,y,z,cellN,cellNF', 12, 9, 'TRI')
g = C.array('x,y,z,cellN,cellNF', 12, 9, 'BAR')
h = C.array('x,y,z,cellN,cellNF', 12, 9, 'NODE')
test.testA([a,b,c,d,e,f,g,h], 2)

# Non structure - array2
a = C.array('x,y,z,cellN,cellNF', 12, 9, 'HEXA', api=2)
#test.testA([a], 4)

test.writeCoverage(100)
