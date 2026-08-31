# basic parametric wing
import Roms.Driver as D
import Geom
import Generator
import Transform

def createPart(name):

    T1 = D.Part(name)

    #======================
    # Create parameters
    #======================
    rootChord = T1.Scalar('rootChord', 1.)

    tipChord = T1.Scalar('tipChord', 1.)

    tipPosx = T1.Scalar('tipPosx', 0.)
    tipPosx.range = [0., 2.]

    tipPosz = T1.Scalar('tipPosz', 2.)
    tipPosz.range = [1.,3.]

    #=================
    # create entities
    #=================

    # Sketch1
    airfoil = Geom.profile("NACA/NACA64-206")

    bbox = Generator.bbox(airfoil)
    grid1 = T1.Grid('grid1', bbox[0:3], bbox[3:], N=(2,2,1))
    spline1 = T1.Spline3('spline1', grid1, mesh=airfoil)
    sketch1 = T1.Sketch('sketch1', [spline1])

    # Sketch2
    airfoil = Geom.profile("NACA/NACA64-206")
    airfoil = Transform.scale(airfoil, 0.5)
    bbox = Generator.bbox(airfoil)
    grid2 = T1.Grid('grid2', bbox[0:3], bbox[3:], N=(2,2,1))
    spline2 = T1.Spline3('spline2', grid2, mesh=airfoil)
    sketch2 = T1.Sketch('sketch2', [spline2])
    sketch2.position.setv(0.,0.,2.)
    T1.Eq(sketch2.position.x, tipPosx)
    T1.Eq(sketch2.position.z, tipPosz)

    # surface
    #surface1 = T1.MergeEdges('surface1', listSketches=[sketch1, sketch2])
    surface1 = T1.Loft('surface1', listSketches=[sketch1, sketch2], h=(1.e-5,0.1,1.e-1))

    #=======================
    # solve and instantiate
    #======================

    # solve
    T1.solve()

    # instantiate
    #T1.instantiate({'tipPosx': 0., 'tipPosz': 2.})

    return T1
