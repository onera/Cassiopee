# - convertPenta2Strand (array) -
import Converter as C
import Generator as G

a = G.cartPenta((0,0,0), (1,1,1), (3,3,3))

b = C.convertPenta2Strand(a)

C.convertArrays2File(b, 'out.msh')
