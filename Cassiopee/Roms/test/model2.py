# - model -
import Roms.Models.nPOD as nPOD
import Converter.PyTree as C

# reload existing model
mod = nPOD.nPOD('density', type='npod')
mod.loadPhi()
mod.loadCoords()
mod.loadPoints()

# instantiate an existing snapshot
t = mod.buildTree(mod.coords[0,:])
C.convertPyTree2File(t, 'out.cgns')

# interpolate in snapshots
t = mod.fetchTree({'Mach':0.75, 'AoA': -1.38})
C.convertPyTree2File(t, 'out.cgns')
