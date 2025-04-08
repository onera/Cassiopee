# - getElementRange (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))

# Get the range of first connectivity
print(Internal.getElementRange(a, number=0))
# Get the range of the first connectivity of TETRA type
print(Internal.getElementRange(a, type='TETRA'))
# Get the range of the connectivity named 'GridElements'
print(Internal.getElementRange(a, name='GridElements'))
