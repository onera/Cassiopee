# driver: parametric circle (derivatives)
import OCC.Driver as D
import Converter

# Create parameter
radius = D.Scalar(1., name='radius')
radius.range = [0.1,10]

# Create circle
circle1 = D.Circle( (0,0,0), radius, name='circle1' )

# Create sketch
sketch1 = D.Sketch([circle1], name='sketch1')

# test
D.DRIVER.solve2()

D.DRIVER.instantiate({'radius': 1.5})

sketch1.writeCAD('out.step')

mesh = sketch1.mesh(0.01, 0.01, 0.01)
D.DRIVER._diff(sketch1, mesh)
Converter.convertArrays2File(mesh, 'out.plt')

#import CPlot, time
#for i in range(50):
#    D.DRIVER.instantiate({'radius': 10-i/10.})
#    mesh = sketch1.mesh(0.01, 0.01, 0.01)
#    CPlot.display(mesh)
#    time.sleep(0.5)