# - convertPyTree2ZoneNames (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
b = G.cartHexa((0.,0.,0.),(0.1,0.1,0.1),(11,11,11)); b[0] = 'cartHexa'
t = C.newPyTree(['Base',a,b])
zones = C.convertPyTree2ZoneNames(t); print(zones)
