# - extCenter2Node (array) -
import Converter as C
import Generator as G

ni = 30; nj = 40; nk = 1
a = G.cart((0,0,0), (1,1,1), (ni,nj,nk))
a = C.initVars(a, '{F}={x}*{x}+{y}*{y}')
ac = C.node2ExtCenter(a)
a = C.extCenter2Node(ac)
C.convertArrays2File([a], "out.plt")
