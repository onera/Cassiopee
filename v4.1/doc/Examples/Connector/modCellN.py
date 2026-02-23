# - modCellN (array) -
import Converter as C
import Generator as G
import Connector as X

a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = C.initVars(a, 'cellN=1.')
a[1][3,22] = 0.
a[1][3,24] = 2.
X._modCellN1(a)
C.convertArrays2File(a, 'out.plt')
