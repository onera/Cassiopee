# - node2Center (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (30,40,1))
a = C.initVars(a, '{ro}=2*{x}+{y}')
ac = C.node2Center(a)
C.convertArrays2File([a,ac], "out.plt")
