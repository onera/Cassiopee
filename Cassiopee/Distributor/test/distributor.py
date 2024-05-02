# - distributor: distribute original blocks amongst computing nodes -

from elsA_user import *
import Converter.Cassiopee as CC
import Converter as C
import Transform as T
import Generator as G
import Distributor 
Distributor.register()

# build blk1
a = G.cart((0.,0.,0.), (0.1,0.1,1.), (10,10,2))
msh = mesh(name='msh')
msh.submit()
CC.convertArray2Mesh(a, msh)
blk = block(name='blk')
blk.set('mesh', 'msh')
blk.submit()

# build blk2
a2 = G.cart((0.2,0.,0.), (0.1,0.1,1.), (10,10,2))
msh2 = mesh(name='msh2')
msh2.submit()
CC.convertArray2Mesh(a2, msh2)
a = T.rotate(a2, (0.2,0.,0.), (0.,0.,1.), 18.)
blk2 = block(name='blk2')
blk2.set('mesh', 'msh2')
blk2.submit()

C.convertArrays2File([a], "mesh.plt", "bin_tp")
C.convertArrays2File([a2], "mesh2.plt", "bin_tp")
