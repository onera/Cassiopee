# - distance2Walls (array) -
import Dist2Walls
import Generator as G
import Converter as C
import Geom as D

# Bloc dont on cherche la distance a la paroi
a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))

# Paroi
sphere = D.sphere((1.2,0.,0.), 0.2, 30)
cellN = C.initVars(sphere,'cellN',1.)
# Calcul de la distance a la paroi
dist = Dist2Walls.distance2Walls(a, [sphere], cellnbodies=[cellN],
                                 loc='centers',type='ortho')
ac = C.node2Center(a)
ac = C.addVars([ac, dist])
C.convertArrays2File([ac], 'out.plt')
