# - getFields (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

z = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
print(C.getFields('GridCoordinates', z))
