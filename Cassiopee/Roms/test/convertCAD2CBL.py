# convertCAD2CBL (array) -
import Roms.Driver as D
import OCC

din = {'sketch1': [1]}

hook = OCC.readCAD('../../OCC/test/out.step')

D.convertCAD2CBL(hook, 'log.py', din)
