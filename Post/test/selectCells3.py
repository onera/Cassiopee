# - selectCells3 (array) -
import Generator as G
import Converter as C
import Post as P

a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = C.node2Center(a)
b = C.initVars(b, '{tag}=0.1*{x}+0.1*{y}')
b = C.extractVars(b, ['tag'])
p = P.selectCells3(a, b)

C.convertArrays2File(p, 'out.plt')
