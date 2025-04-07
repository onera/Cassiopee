# - adapt the bounding box of a point cloud (array) -

import Converter as C
import Generator as G
import Intersector as XOR

a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.convertArray2NGon(a); a = G.close(a)

m = XOR.adaptBox(a, box_ratio=10.)
m = XOR.closeCells(m) # optional : to close the polyhedral cells

C.convertArrays2File([m], 'out.plt')
