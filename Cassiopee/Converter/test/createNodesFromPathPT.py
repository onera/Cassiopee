# - createNodesFromPath (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

t = C.newPyTree(['Base'])
Internal._createNodesFromPath(t, '/CGNSTree/Base/Zone1')
C.convertPyTree2File(t, 'out.cgns')
