# - extractMesh (array) -
import Converter as C
import Post as P
import Generator as G
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
m = G.cart((-5,-5,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
# init by function
m = C.initVars(m, 'F', F, ['x','y','z'])

# Cree un maillage d'extraction
ni2 = 20; nj2 = 20
a = G.cart( (-5.,-5.,0.), (10./(ni2-1),10./(nj2-1),1), (ni2,nj2,2))

# sur une liste - ordre 2 seulement
a2 = P.extractMesh([m], [a], order=2)
test.testA(a2,1)
# Extrait la solution sur le maillage d'extraction
for order in [2,3,5]:
    print('Computing order %d...'%order)
    a2 = P.extractMesh([m], a, order)
    test.testA([a2], order)
