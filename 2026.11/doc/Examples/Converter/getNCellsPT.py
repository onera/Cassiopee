# - getNCells (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,11))
ncells = C.getNCells(a); print(ncells)
#>> 810
