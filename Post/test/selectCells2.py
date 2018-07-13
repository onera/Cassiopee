# - selectCells2 (array) -
import Converter as C
import Generator as G
import Post as P

a = G.cart((0,0,0), (1,1,1), (11,11,11))
tag = C.array('tag', 10, 10, 10); tag = C.initVars(tag, 'tag', 1.)
b = P.selectCells2(a, tag)
C.convertArrays2File([b], 'out.plt')
