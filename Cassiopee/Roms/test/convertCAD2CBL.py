# convertCAD2CBL (array) -
import Roms.Driver as D
import OCC

# create step
hook = OCC.createEmptyCAD()
#OCC._addLine(hook, (0,0,0), (1,0,0))
OCC._addCircle(hook, (0,0,0), (0,0,1), 1.)

# convert step to CBL
din = {'sketch1': [1]}
D.convertCAD2CBL(hook, 'log.py', din)

# test
#import log
