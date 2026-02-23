# - display (array) -
import Generator as G
import CPlot
import time
import Converter as C
import Geom as D

def F(x): return x

# - STRUCTURE 1D -
pts = D.polyline([(0.,0.,0.), (4.,5.,0.), (6.,21.,0.), (18.,2.,0.)])
a = D.bezier(pts, 20)
C._addVars(a, 'p')
a = C.initVars(a, 'p', F, ['x'])

# mesh (mode 0)
CPlot.display([a], mode='mesh', posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display([a], mode='solid'); time.sleep(2)

# scalar field (mode 3)
CPlot.display([a], mode='scalar', scalarField=0); time.sleep(2)

# - BAR -
pts = D.polyline([(0.,0.,0.), (4.,5.,0.), (6.,21.,0.), (18.,2.,0.)])
a = D.bezier(pts, 20)
C._addVars(a, 'p')
a = C.initVars(a, 'p', F, ['x'])
a = C.convertArray2Hexa(a)

# mesh (mode 0)
CPlot.display([a], mode=0, posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display([a], mode=1); time.sleep(2)

# scalar field (mode 3)
CPlot.display([a], mode=3, scalarField=0); time.sleep(2)

# - STRUCTURE 3D -
a = G.cart((0,0,0),(1,1,1),(18,28,3))
a = C.initVars(a, 'p', F, ['x'])

# mesh (mode 0)
CPlot.display([a], mode=0, posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display([a], mode=1); time.sleep(2)

# iso (mode 3)
CPlot.display([a], mode=3, scalarField=0); time.sleep(2)

# - TRI -
a = G.cart((0,0,0),(1,1,1),(18,28,1))
C._addVars(a, 'p')
a = C.initVars(a, 'p', F, ['x'])
a = C.convertArray2Tetra(a)

# mesh (mode 0)
CPlot.display([a], mode=0, posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display([a], mode=1); time.sleep(2)

# scalar field (mode 3)
CPlot.display([a], mode=3, scalarField=0); time.sleep(2)

# - QUAD -
a = G.cart((0,0,0),(1,1,1),(18,28,1))
C._addVars(a, 'p')
a = C.initVars(a, 'p', F, ['x'])
a = C.convertArray2Hexa(a)

# mesh (mode 0)
CPlot.display([a], mode=0, posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display([a], mode=1); time.sleep(2)

# scalar field (mode 3)
CPlot.display([a], mode=3, scalarField=0); time.sleep(2)

# - TETRA -
a = G.cart((0,0,0),(1,1,1),(18,28,3))
C._addVars(a, 'p')
a = C.initVars(a, 'p', F, ['x'])
a = C.convertArray2Tetra(a)

# mesh (mode 0)
CPlot.display([a], mode=0, posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display([a], mode=1); time.sleep(2)

# scalar field (mode 3)
CPlot.display([a], mode=3, scalarField=0); time.sleep(2)

# - HEXA -
a = G.cart((0,0,0),(1,1,1),(18,28,3))
C._addVars(a, 'p')
a = C.initVars(a, 'p', F, ['x'])
a = C.convertArray2Hexa(a)

# mesh (mode 0)
CPlot.display(a, mode=0, posCam=(40,-2,15), posEye=(8,13,0)); time.sleep(2)

# solid (mode 1)
CPlot.display(a, mode=1); time.sleep(2)

# scalar field (mode 3)
CPlot.display(a, mode=3, scalarField=0); time.sleep(2)

# - NODE -
# mesh (mode 0)
a = C.array('x,y,z', 18*28, 0, 'NODE')
for i in range(18):
    for j in range(28):
        C.setValue(a, i+j*18, (i,j,0.))
CPlot.display(a, mode=0); time.sleep(2)
