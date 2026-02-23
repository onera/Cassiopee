# - rotate (PyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.PyTree as C

a = G.cart((0,0,0), (1,1,1), (10,10,2))
# Rotate with an axis and an angle
b = T.rotate(a, (0.,0.,0.), (0.,0.,1.), 30.)
b2 = T.rotate(a, (0.,0.,0.), (0.,0.,1.), 30.,
              vectors=[['centers:MomentumX'], ['centers:MomentumX'], ['centers:MomentumX']])

# Rotate with two axis
c = T.rotate(a, (0.,0.,0.), ((1.,0.,0.),(0,1,0),(0,0,1)),
             ((1,1,0), (1,-1,0), (0,0,1)) )
# Rotate with three angles
c = T.rotate(a, (0.,0.,0.), (0,0,90))

C.convertPyTree2File([a,b,c], 'out.cgns')
