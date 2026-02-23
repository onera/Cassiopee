# - extractPoint : maillage cartesien -
import Converter as C
import Generator as G
import Post as P
#
# Maillage en noeuds
ni = 11; nj = 11; nk = 11;
a = G.cartTetra((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))

# Create a function
def F(x,y,z):
    return 2*x*x*x*x*x + 2.*y*y*z + z*z

(x,y,z) = (0.55,0.38,0.12)
val0 = F(x,y,z)

# init by function
a = C.initVars(a, 'F', F, ['x','y','z'])

cnt = 0
val = P.extractPoint([a], (x,y,z), 3)
err = 0.010897126
if abs(val[0]-val0) > err:
    print('DIFF: reference: '+str(err)+'.')
    print('DIFF: courant: '+str(val[0]-val0)+'.')
