# Test pour Array3
import KCore

import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy

#================================================================
# Array1 - STRUCT
import Generator
a = Generator.cart((0,0,0), (1,1,1), (10,10,10))
#KCore.tester(a)

# Array3 - STRUCT
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = C.getFields('GridCoordinates', a, api=3)[0]
#KCore.tester(a)

#================================================================
# Array1 - MONOBE
a = Generator.cartHexa((0,0,0), (1,1,1), (10,10,10))
#KCore.tester(a)
#print(a[2].shape)
#print(a[2][0,0], a[2][1,0], a[2][2,0])

# Array3 - MONOBE
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
a = C.getFields('GridCoordinates', a, api=3)[0]
#KCore.tester(a)
#print(a)
#print(a[2][0].shape)
#print(a[2][0][0], a[2][0][1], a[2][0][2])

# Array3 - MULTIBE
a1 = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
a2 = G.cartHexa((9,0,0), (1,1,1), (10,10,10))
a = C.mergeConnectivity(a1, a2)
a = C.getFields('GridCoordinates', a, api=3)[0]
KCore.tester(a)
#print(a[2][0].shape)

#=============================================================
# NGON tra la la
