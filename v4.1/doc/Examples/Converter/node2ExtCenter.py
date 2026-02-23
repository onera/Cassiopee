# - node2ExtCenter (array) -
import Converter as C
import Generator as G

ni = 30; nj = 40; nk = 1
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
ac = C.node2ExtCenter(a)
C.convertArrays2File([a,ac], "out.plt")
