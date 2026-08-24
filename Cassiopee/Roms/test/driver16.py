# - driver (array) -
# create a square with shift and sym
import Roms.Driver as D

T1 = D.Part("Part1")

P0 = T1.Point('P0', (-1,-1,0))
P0.x.range = [-2,2]
P0.y.range = [-2,2]
P0.z.range = [0, 1]
P1 = P0.ShiftPoint('P1', (2,0,0))
P2 = P0.SymPoint('P2', 'xz')
P3 = P2.ShiftPoint('P3', (2,0,0))
line1 = T1.Line('line1', P0, P1)
line2 = T1.Line('line2', P1, P3)
line3 = T1.Line('line3', P3, P2)
line4 = T1.Line('line4', P2, P0)

sketch1 = T1.Sketch('sketch1', [line1,line2,line3,line4], h=(0.1,0.1,0.01))

T1.solve()
T1.instantiate({'P0.x': -1.,
                'P0.y': -1.,
                'P0.z': 0.})

sketch1.writeCAD('out.step')

D.DRIVER.plot(T1, sketch1)
