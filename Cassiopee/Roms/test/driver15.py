# driver: parametric cube + volume 3D
import Roms.Driver as D

# Create parameter
hauteur = D.Scalar('hauteur')
hauteur.range = [0.5, 1.5, 0.1]

# Create points
P1 = D.Point('P1', (0,0,0))

P2 = D.Point('P2', (1,0,0))

P3 = D.Point('P3', (1.,1,0))

P4 = D.Point('P4', (0.,1,0))

P5 = D.Point('P5', (0,0,1))
D.Eq(P5.z, hauteur)

P6 = D.Point('P6', (1,0,1))
D.Eq(P6.z, hauteur)

P7 = D.Point('P7', (1.,1,1))
D.Eq(P7.z, hauteur)

P8 = D.Point('P8', (0.,1,1))
D.Eq(P8.z, hauteur)

# Create lines
line1 = D.Line('line1', P1, P2)
line2 = D.Line('line2', P2, P3)
line3 = D.Line('line3', P3, P4)
line4 = D.Line('line4', P4, P1)

line5 = D.Line('line5', P5, P6)
line6 = D.Line('line6', P6, P7)
line7 = D.Line('line7', P7, P8)
line8 = D.Line('line8', P8, P5)

line9 = D.Line('line9', P1, P5)
line10 = D.Line('line10', P2, P6)
line11 = D.Line('line11', P3, P7)
line12 = D.Line('line12', P4, P8)

# Create sketch
sketch1 = D.Sketch('sketch1', [line1, line2, line3, line4], h=[0.01,0.01,0.01])
sketch2 = D.Sketch('sketch2', [line5, line6, line7, line8], h=[0.01,0.01,0.01])
sketch3 = D.Sketch('sketch3', [line1, line10, line5, line9], h=[0.01,0.01,0.01])
sketch4 = D.Sketch('sketch4', [line2, line11, line6, line10], h=[0.01,0.01,0.01])
sketch5 = D.Sketch('sketch5', [line3, line11, line7, line12], h=[0.01,0.01,0.01])
sketch6 = D.Sketch('sketch6', [line4, line12, line8, line9], h=[0.01,0.01,0.01])

# surface
surface1 = D.Fill('surface1', sketch=sketch1, h=[0.1,0.1,0.1])
surface2 = D.Fill('surface2', sketch=sketch2, h=[0.1,0.1,0.1])
surface3 = D.Fill('surface3', sketch=sketch3, h=[0.1,0.1,0.1])
surface4 = D.Fill('surface4', sketch=sketch4, h=[0.1,0.1,0.1])
surface5 = D.Fill('surface5', sketch=sketch5, h=[0.1,0.1,0.1])
surface6 = D.Fill('surface6', sketch=sketch6, h=[0.1,0.1,0.1])
surface = D.Merge('surface', listSurfaces=[surface1,surface2,surface3,surface4,surface5,surface6], h=[0.1,0.1,0.1])

# solve
D.DRIVER.solve()

D.DRIVER.instantiate({'hauteur': 1.5})
#edges = D.MergeEdges(listSketches=[sketch1,sketch2,sketch3,sketch4,sketch5,sketch6])
#edges.writeCAD('out.step')
surface.writeCAD('out.step')
#m = surface1.MeshAsReference()

import CPlot.PyTree as CPlot, time
point = D.DRIVER.walkDOE()
while point is not None:
    D.DRIVER.instantiate(point)
    m = surface.Mesh()
    #m = surface.Dmesh()
    CPlot.display(m)
    point = D.DRIVER.walkDOE()
    time.sleep(0.5)
