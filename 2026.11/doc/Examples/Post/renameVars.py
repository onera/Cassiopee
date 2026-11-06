# - renameVars (array) -
import Converter as C
import Post as P
import Generator as G

ni = 30; nj = 40
m = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,2))
m = C.initVars(m, 'ro', 1.)
m = C.initVars(m, 'rou', 1.)

# Rename a list of variables
m2 = P.renameVars(m, ['ro','rou'], ['Density','MomentumX'])
C.convertArrays2File(m2, "out.plt")
