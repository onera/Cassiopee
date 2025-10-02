# driver: parametric profile
import OCC.Driver as D
import Geom
import Generator
import Converter

# Create parameter
epaisseur = D.Scalar(12., name='epaisseur')
epaisseur.range = [10,15]

# discrete profile
naca = Geom.naca(12, N=51)
bbox = Generator.bbox(naca)

# Create grid
grid1 = D.Grid(bbox[0:3], bbox[3:], N=(3,3,1))
grid1.P[1][2][0].y.range = [0, 5]
D.Eq(epaisseur.s, grid1.P[1][2][0].x.s)

# Create profile
spline1 = D.Spline3( grid1, mesh=naca, name='spline1' )

# Create sketch
sketch1 = D.Sketch([spline1], name='sketch1')

# test
D.DRIVER.solve2()

#grid1.P[1][2][0].y.print()

D.DRIVER.instantiate({'P.1.2.0.y': 0.8})

sketch1.writeCAD('out.step')

mesh = sketch1.mesh(0.01, 0.01, 0.01)
D.DRIVER._diff(sketch1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

D.DRIVER.setDOE({'P.1.2.0.y': 0.1})
D.DRIVER.createDOE('doe.hdf')
D.DRIVER.walkDOE(sketch1, 0.01, 0.01, 0.01)
# reread from doe file
m = D.DRIVER.readSnaphot(0)
Converter.convertArrays2File(m, 'reread.plt')
full = D.DRIVER.readAllSnapshots()
print(full.shape)
#print(full)
D.DRIVER.fullSvd(full)

import sys; sys.exit()

import CPlot, time
for i in range(20):
    D.DRIVER.instantiate({'P.1.2.0.y': 0.3+i/50.})
    mesh = sketch1.mesh(0.01, 0.01, 0.01)
    CPlot.display(mesh)
    time.sleep(0.5)

# build dmesh
mesh1 = D.DRIVER.dmesh(sketch1, mesh, ['P.1.2.0.y'], 0.1)
mesh2 = D.DRIVER.dmesh(sketch1, mesh1, ['P.1.2.0.y'], 0.1)
Converter.convertArrays2File(mesh+mesh1+mesh2, 'out.plt')
