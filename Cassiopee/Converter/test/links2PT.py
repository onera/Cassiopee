# - HDF write with links -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal

a = G.cart((0,0,0),(1,1,1),(50,50,50))
C._initVars(a, 'Density=1.')
t = C.newPyTree(['Base',a])

# Save file with links
links=[['.','coord.hdf','/Base/cart/GridCoordinates','/Base/cart/GridCoordinates']]
C.convertPyTree2File(t, 'main.hdf', links=links)

# Write pointed file
Internal._rmNodeByPath(t, '/Base/cart/FlowSolution')
C.convertPyTree2File(t, 'coord.hdf')
