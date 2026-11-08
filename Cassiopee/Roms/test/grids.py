# test Driver grids
import Roms.Driver as D

T1 = D.Part('test')

v1 = T1.Scalar('v1', 0.)
v1.range = [0., 5., 1.]

grid1 = T1.Grid('grid1', (0,0,0), (1,1,0), (3,3,1))
T1.Eq(grid1.P[0][0][0].y, v1)

grid2 = grid1.Copy('grid2')
grid3 = grid1.Shift('grid3', (1,0,0))
grid4 = grid1.Symm('grid4', plane='yz', center=(0,0,0))
grid5 = grid1.Rotate('grid5', (0,0,0), 30.)

T1.solve()
T1.instantiate({'v1': 1.})

grid1.print()
grid2.print()
grid3.print()
grid4.print()
grid5.print()
