# - extractMesh  sur maillage cylindrique -
import Converter as C
import Post as P
import Generator as G
import Transform as T
import KCore.test as test

# Create a function
def F(x,y,z):
    deg = 1
    if deg == 0 : return 10.
    elif deg == 1 : return x + 2.*y + 3.*z
    elif deg == 2 : return x*x + 2.*y*y + 3*z
    elif deg == 3 : return x*x*y + 2.*y*y*y + 3*z
    elif deg == 4 : return x*x*x*x + 2.*y*y*y*y +z*z
    else : return 2*x*x*x*x*x + 2.*y*y*z + z*z

# Maillage en noeuds
ni = 101; nj = 101; nk = 11
m = G.cylinder((0,0,0), 1., 10.,45., 145., 1., (ni,nj,nk))
m = T.reorder(m, (-1,2,3))
ast = []; ast.append(m)
# init by function
m = C.initVars(m, 'F', F, ['x','y','z'])

# Cree un maillage d'extraction
ni2 = 30; nj2 = 30
a = G.cart( (-1.,2.,0.4), (1./(ni2-1),1./(nj2-1),0.1), (ni2,nj2,2))
ast.append(a)
#
# solution exacte :
c = C.initVars(a, 'F', F, ['x','y','z'])

# sur une liste - ordre 2 seulement
a2 = P.extractMesh([m],[a])
test.testA(a2,1)
# Extrait la solution sur le maillage d'extraction
cnt = 0
for i in [2,3]:
    print('Computing order %d...'%i)
    a2 = P.extractMesh([m], a, i)
    test.testA([a2], i)
