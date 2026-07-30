# driver: parametric cube + volume 3D
import Roms.Driver as D

T1 = D.Part("Part1")

# Create parameter
hauteur = T1.Scalar('hauteur')
hauteur.range = [0.5, 1.5, 0.1]

# Create points
P1 = T1.Point('P1', (0,0,0))

P2 = T1.Point('P2', (1,0,0))

P3 = T1.Point('P3', (1.,1,0))

P4 = T1.Point('P4', (0.,1,0))

P5 = T1.Point('P5', (0,0,1))
T1.Eq(P5.z, hauteur)

P6 = T1.Point('P6', (1,0,1))
T1.Eq(P6.z, hauteur)

P7 = T1.Point('P7', (1.,1,1))
T1.Eq(P7.z, hauteur)

P8 = T1.Point('P8', (0.,1,1))
T1.Eq(P8.z, hauteur)

# Create lines
line1 = T1.Line('line1', P1, P2)
line2 = T1.Line('line2', P2, P3)
line3 = T1.Line('line3', P3, P4)
line4 = T1.Line('line4', P4, P1)

line5 = T1.Line('line5', P5, P6)
line6 = T1.Line('line6', P6, P7)
line7 = T1.Line('line7', P7, P8)
line8 = T1.Line('line8', P8, P5)

line9 = T1.Line('line9', P1, P5)
line10 = T1.Line('line10', P2, P6)
line11 = T1.Line('line11', P3, P7)
line12 = T1.Line('line12', P4, P8)

# Create sketch
sketch1 = T1.Sketch('sketch1', [line1, line2, line3, line4], h=[0.01,0.01,0.01])
sketch2 = T1.Sketch('sketch2', [line5, line6, line7, line8], h=[0.01,0.01,0.01])
sketch3 = T1.Sketch('sketch3', [line1, line10, line5, line9], h=[0.01,0.01,0.01])
sketch4 = T1.Sketch('sketch4', [line2, line11, line6, line10], h=[0.01,0.01,0.01])
sketch5 = T1.Sketch('sketch5', [line3, line11, line7, line12], h=[0.01,0.01,0.01])
sketch6 = T1.Sketch('sketch6', [line4, line12, line8, line9], h=[0.01,0.01,0.01])

# surface
surface1 = T1.Fill('surface1', sketch=sketch1, h=[0.1,0.1,0.1])
surface2 = T1.Fill('surface2', sketch=sketch2, h=[0.1,0.1,0.1])
surface3 = T1.Fill('surface3', sketch=sketch3, h=[0.1,0.1,0.1])
surface4 = T1.Fill('surface4', sketch=sketch4, h=[0.1,0.1,0.1])
surface5 = T1.Fill('surface5', sketch=sketch5, h=[0.1,0.1,0.1])
surface6 = T1.Fill('surface6', sketch=sketch6, h=[0.1,0.1,0.1])
surface = T1.Merge('surface', listSurfaces=[surface1,surface2,surface3,surface4,surface5,surface6], h=[0.1,0.1,0.1])

# solve
T1.solve()

T1.instantiate({'hauteur': 1.5})
#edges = D.MergeEdges(listSketches=[sketch1,sketch2,sketch3,sketch4,sketch5,sketch6])
#edges.writeCAD('out.step')
surface.writeCAD('out.step')
#m = surface1.MeshAsReference()

import CPlot.PyTree as CPlot, time
point = T1.walkDOE()
while point is not None:
    T1.instantiate(point)
    m = surface.Mesh()
    #m = surface.Dmesh()
    CPlot.display(m)
    point = T1.walkDOE()
    time.sleep(0.5)
