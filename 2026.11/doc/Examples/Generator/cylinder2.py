# - cylinder2 (array) -
import Converter as C
import Generator as G
r = G.cart((0.,0.,0.), (0.1, 1., 1.), (11, 1, 1))
teta = G.cart((0.,0.,0.), (0.1, 1., 1.), (11, 1, 1))
z = G.cart((0.,0.,0.), (0.1, 1., 1.), (11, 1, 1))
cyl = G.cylinder2( (0.,0.,0.), 0.5, 1., 360., 0., 10., r, teta, z)
C.convertArrays2File([cyl], "out.plt")
