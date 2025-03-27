# - plot (array) -
import Generator as G
import Transform as T
import Converter as C

# Create a cartesian grid
a = G.cart( (0,0,0), (1,1,1), (10,10,10))

# plot mesh
C.plot( a )

a = T.translate( a, (0.5*1,0,0) )
C.plot( a )

a = T.translate( a, (0.5*1,0,0) )
C.plot( a )
