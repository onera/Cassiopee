# - HDF read/write with links -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Filter as Filter

a = G.cart((0,0,0),(1,1,1),(50,50,50))
t = C.newPyTree(['Base',a])

# Save doesnt follow links
links=[['.','titi.hdf','/Base/cart/GridCoordinates/CoordinateX','/Base/cart/GridCoordinates/CoordinateX']]
C.convertPyTree2File(t, 'tata.hdf', links=links)

# Write links by hand
Internal._rmNodeByPath(t, '/Base/cart/GridCoordinates/CoordinateY')
Internal._rmNodeByPath(t, '/Base/cart/GridCoordinates/CoordinateZ')
C.convertPyTree2File(t, 'titi.hdf')

# full read of tata returning links
LC=[]
t = C.convertFile2PyTree('tata.hdf', links=LC)
print LC

# Read links with skeleton
LC=[]
t = Filter.convertFile2SkeletonTree('tata.hdf', links=LC)
print LC
