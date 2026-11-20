# convertCAD2CBL (array) -
import Roms.Driver as D
import OCC

din = {'sketch1': [1,2,3]}

hook = OCC.readCAD('line.step')

D.convertCAD2CBL(hook, 'log.py', din)

# test
import log
