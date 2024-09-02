# - reorder (array) -
import Generator as G
import Transform as T
import Converter as C

# Structured
a = G.cart((0,0,0),(1,1,1),(5,7,9))
a = T.reorder(a, (3,-2,-1))
C.convertArrays2File(a, 'out1.plt')

# Unstructured
a = G.cartTetra((0,0,0),(1,1,1),(3,3,1))
a = T.reorder(a, (1,))
C.convertArrays2File(a, 'out2.plt')
