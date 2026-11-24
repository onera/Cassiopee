# convertCAD2CBL (array) -
import Roms.Driver as D
import OCC

#din = {'sketch1': [1,2,3,4]}
din = {'loft1': [[1],[4],[5],[6],[2]]}

hook = OCC.readCAD('loft1.step')

D.convertCAD2CBL(hook, 'log.py', din)

# test
import log
