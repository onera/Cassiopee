# - driver (array) -
# create a square with shift and sym
import Roms.Driver as D

P0 = D.Point('P0', (-1,-1,0))
P0.x.range = [-2,2]
P0.y.range = [-2,2]
P0.z.range = [0, 1]
P1 = P0.ShiftPoint('P1', (2,0,0))
P2 = P0.SymPoint('P2', 'yz')
P3 = P2.ShiftPoint('P3', (2,0,0))
line1 = D.Line('line1', P0, P1)
line2 = D.Line('line2', P1, P3)
line3 = D.Line('line3', P3, P2)
line4 = D.Line('line4', P2, P0)

sketch1 = D.Sketch('sketch1', [line1,line2,line3,line4], h=(0.1,0.1,0.01))

D.DRIVER.solve()
D.DRIVER.instantiate({'P0.x': -1.,
                      'P0.y': -1.,
                      'P0.z': 0.})

sketch1.writeCAD('out.step')

D.DRIVER.plot(sketch1)
