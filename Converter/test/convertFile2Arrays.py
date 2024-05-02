# - convertFile2Arrays (arrays) -
import Generator as G
import Converter as C

# Create and save test meshes
cart = G.cart((0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
C.convertArrays2File(cart, 'out.plt')

# Read it
A = C.convertFile2Arrays('out.plt'); print(A)
#>> [['x,y,z', array([[ 0. ,  0.1,  0.2,  ...]]), 11,11,2]]
