# - rmVars (array) -
import Converter as C
import Generator as G
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = C.addVars(a, 'Density')
b = C.addVars(a, 'Alpha')
b = C.rmVars(b, 'Density')
C.convertArrays2File(b, 'out.plt')
