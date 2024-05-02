# - addVars (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (10,10,11))
b = C.array('cell', a[2], a[3], a[4])
c = C.array('t,u', a[2], a[3], a[4])
f = C.addVars([a, b, c])
C.convertArrays2File(f, 'out.plt')
