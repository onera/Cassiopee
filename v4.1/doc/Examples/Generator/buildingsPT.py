# - Buildings (pyTree) -
import Generator.Buildings as Buildings
import Converter.PyTree as C

# les dimensions sont en metres
#base = Buildings.square((0.,0.,0.), (100.,50.,0.))
#domain = Buildings.domain(base, height=20., h=1., hp=0.01)
base = Buildings.square((40.,10.,0.), (60.,20.,0.))
building1 = Buildings.building(base, height=10., h=1., hp=0.01)
#t = C.newPyTree(['Fond',domain,'Building1',building1])
t = C.newPyTree(['Building1',building1])
C.convertPyTree2File(t, 'out.cgns')