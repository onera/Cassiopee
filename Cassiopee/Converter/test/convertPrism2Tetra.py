# - convertArray2Tetra (array) -
import Converter as C
import Generator as G

a = G.cartPenta((0.,0.,0.), (1,1,1), (10,10,3))
a = C.convertArray2Tetra(a)
C.convertArrays2File([a], 'out.plt', 'bin_tp')
