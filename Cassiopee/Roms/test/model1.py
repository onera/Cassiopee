# - model -
import Roms.DB.DataBase as DataBase
import Roms.Models.POD as POD

db = DataBase.DataBase('NACA1.db', mode='r')
q = db.query()
A = db.fetchMatrix(q, variables=['centers:Density'])

mod = POD.POD('density', type='pod')
mod.buildPhi(A)
mod.savePhi()
mod.buildAndSaveCoeffs(A)