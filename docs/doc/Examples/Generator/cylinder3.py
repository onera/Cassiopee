# - cylinder3 (array) -
import Generator as G
import Converter as C
teta = G.cart((0.,0.,0.), (0.1, 1., 1.), (11, 1, 1))
xz = G.cart((0.1,0.,0.), (0.1,1.,0.2), (20, 1, 30))
cyl = G.cylinder3( xz, 0., 90., teta)
C.convertArrays2File([cyl], 'out.plt')
