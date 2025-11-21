# dataBase
import Roms.DB.DataBase as DataBase

# create dataBase
db = DataBase.DataBase("NACA1", parameters=['Mach', 'AoA'])

# register without zone
db.register(1, "sample1", {'Mach':0.5, 'AoA':0.0}, variables=["Pressure"], data=None)
db.register(2, "sample2", {'Mach':0.6, 'AoA':0.0}, variables=["Pressure"], data=None)
db.register(3, "sample3", {'Mach':0.7, 'AoA':1.0}, variables=["Pressure"], data=None)

# Create data (pytree)
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
z = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(z, '{Pressure}={CoordinateX}')
t = C.newPyTree(['Base', z])
db.registerReference(t, "ref1")

# register data
db.register(1, "sample1", {'Mach':0.5, 'AoA':0.0}, "ref1", data=t)
db.register(2, "sample2", {'Mach':0.6, 'AoA':0.0}, "ref1", data=t)
db.register(3, "sample3", {'Mach':0.7, 'AoA':1.0}, "ref1", data=t)


# catalog
q = db.query()
print(q)

# query
q = db.query("id = 1")
print(q)
q = db.query("Mach = 0.7 AND AoA >= 0.0")
print(q)
ts = db.fetchTree(q, mode=0)
Internal.printTree(ts)

param = db.fetchParam(q)
print(param)

matrix = db.fetchMatrix(q, variable='Pressure')
print(matrix)