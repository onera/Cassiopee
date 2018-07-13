# - rotate (PyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,2))
# Rotate with an axis and an angle
b = T.rotate(a, (0.,0.,0.), (0.,0.,1.), 30.); b[0] = 'cartRot1'
# Rotate with two axis
c = T.rotate(a, (0.,0.,0.), ((1.,0.,0.),(0,1,0),(0,0,1)),
             ((1,1,0), (1,-1,0), (0,0,1)) ); c[0] = 'cartRot2'
# Rotate with three angles
c = T.rotate(a, (0.,0.,0.), (0,0,90)); c[0] = 'cartRot3'
C.convertPyTree2File([a,b,c], 'out.cgns')
