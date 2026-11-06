# - dLength (array) -
import Geom as D
import Converter as C

a = D.line((0,0,0), (1,0,0), N=2)
print(D.getLength(a))
a = C.addVars(a, ['dx','dy','dz'])

b = D.dLength(a)
print(b)
