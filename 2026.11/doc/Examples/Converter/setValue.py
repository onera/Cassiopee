# - setValue (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1,1,1), (5,5,1))

# Set point (1,1,1) with value x=0.1, y =0.1, z=1.
C.setValue(a, (1,1,1), [0.1,0.1,1.]); print(a)

# Same thing with a global index
C.setValue(a, 0, [0.1,0.1,1.])
