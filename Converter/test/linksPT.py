# - HDF read/write with links -

import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter

a = G.cart((0,0,0),(1,1,1),(50,50,50))
t = C.newPyTree(['Base',a])
C.convertPyTree2File(t, 'coord.hdf')
C._initVars(t, 'Density=1.')

# Save file with links
links=[['.','coord.hdf','/Base/cart/GridCoordinates','/Base/cart/GridCoordinates']]
C.convertPyTree2File(t, 'main.hdf', links=links)

# full read of main returning links
LC=[]
t = C.convertFile2PyTree('main.hdf', links=LC); print(LC)
#>> [['.', './coord.hdf', '/Base/cart/GridCoordinates', '/Base/cart/GridCoordinates']]

# Read links with skeleton
LC=[]
t = Filter.convertFile2SkeletonTree('main.hdf', links=LC); print(LC)
#>> [['.', './coord.hdf', '/Base/cart/GridCoordinates', '/Base/cart/GridCoordinates']]
