# - rotate (array) -
import Generator as G
import Transform as T
import Converter as C

a = G.cart((0,0,0), (1,1,1), (10,10,1))
# Rotate with an axis and an angle
b = T.rotate(a, (0.,0.,0.), (0.,0.,1.), 30.)
# Rotate with axis transformations
c = T.rotate(a, (0.,0.,0.), ((1.,0.,0.),(0,1,0),(0,0,1)),
             ((1,1,0), (1,-1,0), (0,0,1)) )
# Rotate with three angles
d = T.rotate(a, (0.,0.,0.), (90.,0.,0.))
C.convertArrays2File([a,d], 'out.plt')
