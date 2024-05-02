# - getZoneDim (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
dim = Internal.getZoneDim(a); print(dim)
#>> ['Structured', 10, 10, 10, 3]

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
dim = Internal.getZoneDim(a); print(dim)
#>> ['Unstructured', 1000, 3645, 'TETRA', 3]
