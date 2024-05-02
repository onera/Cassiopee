# - getArgMin (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (10,10,10))
a = C.initVars(a, '{F}=({x}-3)*({x}-3)+{y}+{z}')
argmin = C.getArgMin(a, 'F'); print(argmin)
# Corresponds to x, y, z, F for the node where F is mininum
#>> [0.0, 0.0, 0.0, 0.0]
