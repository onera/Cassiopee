import Generator as G
import Converter as C

a = G.cart( (0,0,0), (1,1,1), (5,5,5))
a = C.initVars(a, '{F} = logical_and({x}>1,{y}>1)')
C.convertArrays2File([a], 'out.plt')
