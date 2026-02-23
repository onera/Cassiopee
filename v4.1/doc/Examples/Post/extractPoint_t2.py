# - extractPoint (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test
#
# Maillage en noeuds
ni = 11; nj = 11; nk = 1;
a = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1.), (ni,nj,nk))

# Create a function
def F(x,y,z):
    return 2*x + 2.*y + z

ref = [F(0.55,0.38,0.)]
# init by function
a = C.initVars(a, 'F', F, ['x','y','z'])
val = P.extractPoint([a], (0.55, 0.38, 0.))
# Une fonction lineaire doit etre interpolee exactement
print("Test1... done.")
for i in range(len(val)):
    if abs(val[i]-ref[i]) > 1.e-10:
        print('DIFF: reference: '+str(ref[i])+'.')
        print('DIFF: courant: '+str(val[i])+'.')
