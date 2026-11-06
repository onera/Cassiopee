# - getIndexField (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,1,1))
n = C.getIndexField(a)
a = C.addVars([a, n])
C.convertArrays2File(a, 'out.plt')
