# - randomizeVar (array) -
import Converter as C
import Generator as G
a = G.cart((0,0,0),(1,1,1),(11,11,1))
a = C.initVars(a, '{F}={x}*{y}')
b = C.randomizeVar(a,'F',0.1,0.5)
C.convertArrays2File([a,b],"out.plt")
