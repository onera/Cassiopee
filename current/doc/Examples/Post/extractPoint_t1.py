# - extractPoint (array) -
import Converter as C
import Generator as G
import Post as P
import KCore.test as test

# Tests sur maillage cartesien
ni = 10; nj = 10; nk = 10;
a = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))

# Create a function
def F(x,y,z):
    if deg == 0 : return 10.
    elif deg == 1: return x + 2.*y + 3.*z
    elif deg == 2: return x*x + 2.*y*y + 3*z
    elif deg == 3: return x*x*y + 2.*y*y*y + 3*z
    elif deg == 4: return x*x*x*x + 2.*y + z*z
    elif deg == 5: return 2*x*x*x*x*x + 2.*y*y*z + z*z
    else:
        print('Error: unknown order.')
        import sys; sys.exit()

deg = 3
a = C.initVars(a, 'F', F, ['x','y','z'])

# Interpole 1 point a l'ordre 2
val = P.extractPoint([a], (0.55, 0.38, 0.12), 2)
ref = [0.5918312757201647]
print("Test1... done.")
for i in range(len(val)):
    if abs(val[i]-ref[i]) > 1.e-10:
        print('DIFF: reference: '+str(ref[i])+'.')
        print('DIFF: courant: '+str(val[i])+'.')

# Interpole 1 point a l'ordre 3
val = P.extractPoint([a], (0.55, 0.38, 0.12), 3)
ref = [0.58564300411522652]
print("Test2... done.")
for i in range(len(val)):
    if abs(val[i]-ref[i]) > 1.e-10:
        print('DIFF: reference: '+str(ref[i])+'.')
        print('DIFF: courant: '+str(val[i])+'.')

# Interpole 1 point a l'ordre 5
val = P.extractPoint([a], (0.55, 0.38, 0.12), 5)
ref = [0.58469400000224414]
print("Test3... done.")
for i in range(len(val)):
    if abs(val[i]-ref[i]) > 1.e-10:
        print('DIFF: reference: '+str(ref[i])+'.')
        print('DIFF: courant: '+str(val[i])+'.')

# En dehors
val = P.extractPoint([a], (2,0,0), 2)

# Interpole une liste de points
val = P.extractPoint([a], [(0.1,0.1,0.1), (0.2,0.2,0.3)], 2)

# Maillage en noeuds
ni = 11; nj = 11; nk = 11;
a = G.cart((0,0,0), (1./(ni-1),1./(nj-1),1./(nk-1)), (ni,nj,nk))

val0 = F(0.55,0.38,0.12)

# init by function
C._addVars(a, 'F')
a = C.initVars(a, 'F', F, ['x','y','z'])

cnt = 0
err5 = [0.009497125, 0.001978125, 2.81249995571e-05]
for i in [2,3,5]:
    val = P.extractPoint([a], (0.55, 0.38, 0.12), i)
    print("Test order "+str(i)+"... done.")
    if abs(val[0]-val0) > err5[cnt]:
        print('DIFF: reference: '+str(err5[cnt])+'.')
        print('DIFF: courant: '+str(val[0]-val0)+'.')
    cnt += 1

cnt = 0
hook = C.createHook([a], function='extractMesh')
for i in [2,3,5]:
    val = P.extractPoint([a], (0.55, 0.38, 0.12), i)
    print("Test order "+str(i)+"... done.")
    if abs(val[0]-val0) > err5[cnt]:
        print('DIFF: reference: '+str(err5[cnt])+'.')
        print('DIFF: courant: '+str(val[0]-val0)+'.')
    cnt += 1
C.freeHook(hook)
