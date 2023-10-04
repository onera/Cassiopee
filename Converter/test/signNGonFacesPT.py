# - signNGonFaces (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
ar = C.getFields('GridCoordinates', a, api=3)[0]

import Converter
Converter.converter.signNGonFaces(ar)
