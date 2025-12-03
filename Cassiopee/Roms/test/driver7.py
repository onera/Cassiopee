# driver: parametric revolve surface
import Roms.Driver as D
import Converter

# Create parameter
epaisseur = D.Scalar('epaisseur', 1.)
epaisseur.range = [0.5,2.]

# create Points
P1 = D.Point('P1', (0.1,0,0))
P2 = D.Point('P2', (1.,0.,1))
P2.x.range = [0.5,2.]
P3 = D.Point('P3', (0.1,0.,2))

D.Eq(epaisseur.s, P2.x.s)

# Create profile
spline1 = D.Spline1('spline1', [P1,P2,P3])

# Create sketch 1
sketch1 = D.Sketch('sketch1', [spline1])

# surface
surface1 = D.Revolve('surface1', sketch1, center=(0,0,0), axis=(0,0,1), angle=90.)

# test
D.DRIVER.solve2()

D.DRIVER.instantiate({'epaisseur': 1.5})

surface1.writeCAD('out.step')

mesh = surface1.mesh(0.05, 0.05, 0.1)
D.DRIVER._diff(surface1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

#import CPlot, time
#for i in range(50):
#    D.DRIVER.instantiate({'P2.x': 0.5+i/50.})
#    mesh = surface1.mesh(0.01, 0.01, 0.01)
#    CPlot.display(mesh)
#    time.sleep(0.5)