# dataBase
import Roms.DB.DataBase as DataBase

# create dataBase
db = DataBase.DataBase("NACA1", parameters=['Mach', 'AoA'])
q = db.query()
print(q)

# register without zone
db.register("sample", {'Mach':0.5, 'AoA':0.0}, variables=["Pressure"], data=None)
db.register("sample", {'Mach':0.6, 'AoA':0.0}, variables=["Pressure"], data=None)
db.register("sample", {'Mach':0.7, 'AoA':1.0}, variables=["Pressure"], data=None)

# Create data (pytree)
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
z = G.cart((0,0,0), (1,1,1), (10,10,10))
C._initVars(z, '{Pressure}={CoordinateX}')
t = C.newPyTree(['Base', z])
db.registerReference(t, "ref1")

# register data (replace)
db.register("sample", {'Mach':0.5, 'AoA':0.0}, "ref1", data=t)
db.register("sample", {'Mach':0.6, 'AoA':0.0}, "ref1", data=t)
db.register("sample", {'Mach':0.7, 'AoA':1.0}, "ref1", data=t)

# query catalog
print(db.query())

# query items
q = db.query("id = 1")
print(q)
q = db.query("Mach = 0.7 AND AoA >= 0.0")
print(q)
ts = db.fetchTree(q, mode=0)[0]
#C.convertPyTree2File(ts, 'out.cgns')
Internal.printTree(ts)

# fetch param vector
param = db.fetchParam(q)
print(param)

# fetch param points
points = db.fetchPoint(q)
print(db.exist(points[0]))
print(points)

# fetch matrix
matrix = db.fetchMatrix(q, variable='Pressure')
print(matrix)

# delete query
q = db.query("Mach = 0.7 AND AoA = 1.0")
db.delete(q)
