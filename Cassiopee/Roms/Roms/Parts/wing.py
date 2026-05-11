# parametric wing
import Roms.Driver as D
import Geom
import Generator

#======================
# Create parameters
#======================

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
bbox = Generator.bbox(airfoil)
grid2 = D.Grid('grid2', bbox[0:3], bbox[3:], N=(2,2,1))
spline2 = D.Spline3('spline2', grid2, mesh=airfoil)
sketch2 = D.Sketch('sketch2', [spline2])
sketch2.position.setv(0.2,0.,2.)

# surface
#surface1 = D.MergeEdges('surface1', listSketches=[sketch1, sketch2])
surface1 = D.Loft('surface1', listSketches=[sketch1, sketch2])

#=======================
# solve and instantiate
#======================

# solve
D.DRIVER.solve()

# instantiate
D.DRIVER.instantiate({})

# export
surface1.writeCAD('out.step')
