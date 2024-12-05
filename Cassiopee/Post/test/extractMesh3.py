# - extractMesh (array) -
import Converter as C
import Post as P
import Generator as G
import Transform as T

# Create a function
def F(x,y,z):
    deg = 1
    if deg == 0 :
        return 10.
    elif deg == 1 :
        return x + 2.*y + 3.*z
    elif deg == 2 :
        return x*x + 2.*y*y + 3*z
    elif deg == 3 :
        return x*x*y + 2.*y*y*y + 3*z
    elif deg == 4 :
        return x*x*x*x + 2.*y*y*y*y +z*z
    elif deg == 5 :
        return 2*x*x*x*x*x + 2.*y*y*z + z*z
    else:
        print('Error : unknown degree of polynomials')
        import sys
        sys.exit()

# Maillage en noeuds
ni = 101; nj = 101; nk = 11
m = G.cylinder((0,0,0), 1., 10.,45., 145., 1., (ni,nj,nk))
m = T.reorder(m, (-1,2,3))

ars = []
ars.append(m)
# init by function
m = C.addVars(m, 'F')
m = C.initVars(m, 'F', F, ['x','y','z'])

# Cree un maillage d'extraction
ni2 = 30; nj2 = 30
a = G.cart( (-1.,2.,0.4), (1./(ni2-1),1./(nj2-1),0.1), (ni2,nj2,2))
ars.append(a)
C.convertArrays2File(ars, "out.plt")

# Extrait la solution sur le maillage d'extraction
ordret = [2,3,5]
for ordre in ordret:
    a2 = P.extractMesh([m], a, ordre)
    # solution exacte :
    c = C.addVars(a, 'F')
    c = C.initVars(c, 'F', F, ['x','y','z'])
    #
    e = C.diffArrays([a2], [c])
    print('ordre %d=%f'%(ordre, C.normL0(e[0],"DF")))

# ordre  2  =  2.6645352591e-15
# ordre  3  =  5.09975549789e-06
# ordre  5  =  5.88661031031e-09
