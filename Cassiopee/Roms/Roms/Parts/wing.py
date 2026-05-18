# basic parametric wing
import Roms.Driver as D
import Geom
import Generator
import Transform
import Converter

#======================
# Create parameters
#======================
rootChord = D.Scalar('rootChord', 1.)

tipChord = D.Scalar('tipChord', 1.)

tipPosx = D.Scalar('tipPosx', 0.)
tipPosx.range = [0., 2.]

tipPosz = D.Scalar('tipPosz', 2.)
tipPosz.range = [1.,3.]

#=================
# create entities
#=================

# Sketch1
airfoil = Geom.profile("NACA/NACA64-206")

bbox = Generator.bbox(airfoil)
grid1 = D.Grid('grid1', bbox[0:3], bbox[3:], N=(2,2,1))
spline1 = D.Spline3('spline1', grid1, mesh=airfoil)
sketch1 = D.Sketch('sketch1', [spline1])

# Sketch2
airfoil = Geom.profile("NACA/NACA64-206")
airfoil = Transform.scale(airfoil, 0.5)
bbox = Generator.bbox(airfoil)
grid2 = D.Grid('grid2', bbox[0:3], bbox[3:], N=(2,2,1))
spline2 = D.Spline3('spline2', grid2, mesh=airfoil)
sketch2 = D.Sketch('sketch2', [spline2])
sketch2.position.setv(0.,0.,2.)
D.Eq(sketch2.position.x, tipPosx)
D.Eq(sketch2.position.z, tipPosz)

# surface
#surface1 = D.MergeEdges('surface1', listSketches=[sketch1, sketch2])
surface1 = D.Loft('surface1', listSketches=[sketch1, sketch2], h=(1.e-5,0.1,1.e-1))

#=======================
# solve and instantiate
#======================

# solve
D.DRIVER.solve()

# instantiate
D.DRIVER.instantiate({'tipPosx': 0., 'tipPosz': 2.})

# export CAD and mesh
surface1.writeCAD('out.step')

m = surface1.mesh(method=0)
Converter.convertArrays2File(m, 'out.plt')

D.DRIVER.plot(surface1)
