# - isNamePresent (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(a, 'F', 1.)

b = G.cart((0,0,0), (1,1,1), (50,50,50))
C._initVars(b, 'centers:G', 2.)

print(C.getVarNames([a, b]))
#>> [['x', 'y', 'z', 'F'], ['x', 'y', 'z', 'G']]

print(C.isNamePresent(a, 'F'))
#>> 1
print(C.isNamePresent([a, b], 'F'))
#>> 0
print(C.isNamePresent([a, b], 'K'))
#>> -1
