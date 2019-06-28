# - adapts a cells with respect to b points (array) -
import Intersector as XOR
import Converter as C
import Generator as G

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)
#C.convertArrays2File([a], 'a.plt')

nodal_vals = C.array('cellN',4,4,4)
nodal_vals = C.initVars(nodal_vals, 'cellN', 2.)
nodal_vals = C.center2Node(nodal_vals)

m = XOR.adaptCellsNodal(a, nodal_vals)

m = XOR.closeOctalCells(m[0])
C.convertArrays2File([m], 'out.plt')


