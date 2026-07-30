# new driver using parts
import Roms.Driver as D

T1 = D.DRIVER.Part('square1', type='Square')
T1.instantiate({'width':1., 'height':2.})

T2 = D.DRIVER.Part('circle1', type='Circle')
T2.instantiate({'radius':1.})

#D.DRIVER.print()
D.DRIVER.writeCAD('out.step')
D.DRIVER.plot()
