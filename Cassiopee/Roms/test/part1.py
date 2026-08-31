# new driver using parts
import Roms.Driver as D

T1 = D.DRIVER.Part('square1', type='Square')
T1.instantiate({'width':1., 'height':2.})

#T2 = D.DRIVER.Part('circle1', type='Circle')
#T2.instantiate({'radius':1.})

#T1 = D.DRIVER.Part('nozzle1', type='nozzle2D')
#T1.instantiate({'P1.x': 0.,
#                'Hi': 1.,
#                'Ht': 0.2,
#                'He': 1.,
#                'Lc': 2.,
#                'Ld': 2.})

#T1 = D.DRIVER.Part('nozzle1', type='nozzle3D')
#T1.instantiate({'P1.x': 0.,
#                'Hi': 1.,
#                'Ht': 0.2,
#                'He': 1.,
#                'Lc': 2.,
#                'Ld': 2.})

#T1 = D.DRIVER.Part('wedge1', type='wedge')
#T1.instantiate({'Lf': 2,
#                'Lb': 0.1,
#                'Ht': 0.1,
#                'Hb': 0.1,
#                'Ls': 1.})

#T1 = D.DRIVER.Part('wing1', type='wing')
#T1.instantiate({'tipPosx': 0., 'tipPosz': 2.})

#T1 = D.DRIVER.Part('wing1', type='missile')
#T1.instantiate({'radius1': 10., 'radius2': 4.})

#T1 = D.DRIVER.Part('wing1', type='ImbricatedSquares')
#T1.instantiate({'P0#2.y':1., 'P0#2.x':1., 'depth':1., 'angle':20., 'width':1.})
#D.DRIVER.print()

D.DRIVER.writeCAD('out.step')
D.DRIVER.plot()
