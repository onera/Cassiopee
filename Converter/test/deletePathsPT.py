# - deletePaths (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Filter as Filter

t = C.newPyTree(['Base'])
C.convertPyTree2File(t, 'out.adf')

# Delete paths
Filter.deletePaths('out.adf', 'CGNSTree/Base')
