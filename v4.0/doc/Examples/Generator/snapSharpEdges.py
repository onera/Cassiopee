# - snapSharpEdges (array) -
import Generator as G
import Converter as C
import Post as P
import Geom as D

# Enforce polyline define by s in b
s = D.polyline([(0.2,0,0),(1,1,0),(2.5,1,0),(0.2,0,0)])
s = C.initVars(s, 'indic', 0)
h = 0.1
ni = 30; nj = 20; nk=1
b = G.cartHexa((-0.5, -0.5, 0), (h, h, 1.), (ni,nj,nk))
b = C.initVars(b, 'indic', 0)
b = G.snapSharpEdges(b, [s], h)
c = C.converter.convertQuad2Tri(b)
C.convertArrays2File([b,c, s], 'out.plt')

# Same with smooth
#c = T.smooth(c, eps=0.5, niter=5,
#             fixedConstraints=s)
#            #projConstraints=s)C

# Enforce all constraints (must be over-refined)
s = D.circle((0,0,0), R=1, N=400)
s = C.initVars(s, 'indic', 0)
h = 0.3
ni = 10; nj = 20; nk=1
b = G.cartHexa((-1.5, -1.5, 0), (h, h, 1.), (ni,nj,nk))
b = C.initVars(b, 'indic', 0)
b = G.snapSharpEdges(b, [s], 0.1*h)
c = C.converter.convertQuad2Tri(b)
C.convertArrays2File([b,c,s], 'out.plt')
import sys; sys.exit()

# Idem external
h = 0.3
ni = 10; nj = 10; nk=1
s = G.cartHexa((-1.6,-1.6,0), (h/2,h/2,1.), (2*ni+2,2*nj+2,nk))
s = P.exteriorFaces(s)
s = C.initVars(s, 'indic', 0)
b = G.cartHexa((-1.5, -1.5, 0), (h, h, 1.), (ni,nj,nk))
b = C.initVars(b, 'indic', 0)
b = G.snapSharpEdges(b, [s], h*0.1)
C.convertArrays2File([b,s], 'out.plt')
