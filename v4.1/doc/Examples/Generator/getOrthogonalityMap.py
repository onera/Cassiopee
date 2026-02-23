# - getOrthogonalityMap (array) -
import Generator as G
import Converter as C

a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
ac = C.node2Center(a)
ortho = G.getOrthogonalityMap(a)
ortho = C.addVars([ac,  ortho])
C.convertArrays2File([ortho], "out.plt")
