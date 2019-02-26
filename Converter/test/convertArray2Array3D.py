# - convertArrays2Arrays3D -
import Generator as G
import Converter.Array3D

a = G.cart((0,0,0), (0.1, 0.2, 1.), (11, 4, 1))
b = Converter.Array3D.convertArrays2Arrays3D([a]); print(b)
#>> [[['x', 'y', 'z'], [array([[[ 0. ], ...]]]
