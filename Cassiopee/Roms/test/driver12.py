# driver: parametric sketch for naca4 digits
# + volume meshing
import Roms.Driver as D

T1 = D.Part("Part1")

# Thickness
xx = T1.Scalar('xx')
xx.range = [10., 14., 1.]

# max camber
M = T1.Scalar('M')
M.range = [0., 5., 1.]

# camber position
P = T1.Scalar('P')
P.range = [0., 5., 1.]

# Contrainte
T1.Le(M, P)

# Create points for exterior domain
P1 = T1.Point('P1', (-20,-20,0))

P2 = T1.Point('P2', (-20,20,0))

P3 = T1.Point('P3', (20,20,0))

P4 = T1.Point('P4', (20,-20,0))

# Create lines
line1 = T1.Line('line1', P1, P2)
line2 = T1.Line('line2', P2, P3)
line3 = T1.Line('line3', P3, P4)
line4 = T1.Line('line4', P4, P1)

# create profile
naca1 = T1.Naca('naca1', M, P, xx)

# Create sketch
sketch1 = T1.Sketch('sketch1', [naca1], h=[0.005,0.005,0.005])

sketch2 = T1.Sketch('sketch2', [line1,line2,line3,line4], h=[3.,3.,3.])

# Create volume
vol1 = D.Volume2D('vol1', [sketch1, sketch2], orders=[+1,-1])

# solve
T1.solve()

T1.instantiate({'M':0., 'P': 0., 'xx': 12.})
m = vol1.MeshAsReference()

import CPlot.PyTree as CPlot
import time
pt = T1.walkDOE()
while pt is not None:
    #m = vol1.Mesh()
    m = vol1.Dmesh()
    CPlot.display(m)
    CPlot.setState(message="M=%d P=%d xx=%d"%(M.v,P.v,xx.v))
    pt = T1.walkDOE()
    time.sleep(0.5)
