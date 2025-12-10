# - model -
import Roms.DB.DataBase as DataBase
import Roms.Models.nPOD as nPOD
import Converter.PyTree as C

db = DataBase.DataBase('NACA1.db', mode='r')
q = db.query()
A = db.fetchMatrix(q, variables=['centers:Density'])
W = db.fetchW("ref1")
P = db.fetchPointVector(q)

mod = nPOD.nPOD('density', type='npod')
mod.set(A, W, P, db=db, ref="ref1", parameters=db.parameters, variables=['centers:Density'])

# Build phi
mod.buildPhi()
mod.savePhi()
mod.buildCoords()
mod.saveCoords()
mod.savePoints()

# instantiate an existing snapshot
t = mod.buildTree(mod.coords[0,:])
C.convertPyTree2File(t, 'out.cgns')

# interpolate in snapshot
t = mod.fetchTree({'Mach':0.75, 'AoA': -1.38})
C.convertPyTree2File(t, 'out.cgns')
