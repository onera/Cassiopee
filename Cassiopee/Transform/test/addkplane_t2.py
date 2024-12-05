# - addkplane w/ Flowfield values -
import Transform.PyTree as T
import KCore.test as test
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.), (0.1,0.1,1), (10,10,2))
t = C.newPyTree(['CARTESIAN', a])
C._initVars(t,'{RandomVal}={CoordinateY}')
t = C.node2Center(t, 'RandomVal')
T._addkplane(t, 9)
T._contract(t, (0,0,0), (1,0,0), (0,1,0), 0.1)

test.testT(t ,1)

#C.convertPyTree2File(t,'t_check.cgns')
