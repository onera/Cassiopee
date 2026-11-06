# - renameVars (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.addVars(m, ['Density', 'centers:MomentumX'])

# Rename a list of variables
m2 = P.renameVars(m, ['Density', 'centers:MomentumX'], ['Density_M', 'centers:MomentumX_M'])

C.convertPyTree2File(m2, 'out.cgns')
