# - getValue of a node (pyTree) -
import Generator.PyTree as G
import Converter.Internal as Internal

# Structured array
a = G.cart((0,0,0), (1., 0.5,1.), (40,50,20))

# Get value stored in a zone node
print(Internal.getValue(a))
#>> [[40 39  0] [50 49  0] [20 19  0]]

# Get type of a zone (from ZoneType node)
node = Internal.getNodeFromName(a, 'ZoneType')

# Print node[1], which is a numpy array
print(node[1])
#>> array(['S', 't', 'r', 'u', 'c', 't', 'u', 'r', 'e', 'd']

# getValue, return a string in this case
print(Internal.getValue(node))
#>> Structured
