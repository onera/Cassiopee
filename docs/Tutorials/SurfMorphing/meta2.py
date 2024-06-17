# Linear transformation of two surfaces
import Converter as C
import Dist2Walls
import CPlot
import Generator as G
import Post as P

a = C.convertFile2Arrays('case2.plt')
shape1 = a[0]
shape2 = a[1]

density = 0.1

# Grid autour de a
BB = G.bbox(a)
xmin = BB[0]; ymin = BB[1]; zmin = BB[2]
xmax = BB[3]; ymax = BB[4]; zmax = BB[5]
ni = density*(xmax-xmin); nj = density*(ymax-ymin);
nk = density*(zmax-zmin)
if (ni < 2): ni = 2
if (nj < 2): nj = 2
if (nk < 2): nk = 2

hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1) ; hk = (zmax-zmin)/(nk-1)
h = min(hi, hj) ; h = min(h, hk)

ni = int((xmax-xmin)/h)+7 ; nj = int((ymax-ymin)/h)+7
nk = int((zmax-zmin)/h)+7

b = G.cart( (xmin-3*h, ymin-3*h, zmin-3*h), (h, h, h), (ni,nj,nk) )

# Distance initiale (shape1)
d1 = Dist2Walls.distance2Walls(b, [shape1], type='ortho', loc='nodes', signed=1)
d1[0] = 'Dist1'

# Distance finale (shape2)
d2 = Dist2Walls.distance2Walls(b, [shape2], type='ortho', loc='nodes', signed=1)
d2[0] = 'Dist2'

b = C.addVars([b, d1])
b = C.addVars([b, d2])

# Lets groove
val = 0.*h

for j in range(10):
    
    for i in range(200):
        alpha = i/199.
        b = C.initVars(b, 'Dist = '+str(alpha)+'*{Dist1}+'+str((1.-alpha))+'*{Dist2}')
        iso = P.isoSurfMC([b], 'Dist', value=val)
        CPlot.display(iso, mode=0, meshStyle=0)
        #CPlot.setState(stereo=1)

    for i in range(200):
        alpha = i/199.
        b = C.initVars(b, 'Dist = '+str(1.-alpha)+'*{Dist1}+'+str(alpha)+'*{Dist2}')
        iso = P.isoSurfMC([b], 'Dist', value=val)
        CPlot.display(iso, mode=0, meshStyle=0)
