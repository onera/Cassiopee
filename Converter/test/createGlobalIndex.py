# - createGlobalIndex (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
C._createGlobalIndex(a)
C.convertArrays2File(a, 'out.plt')
