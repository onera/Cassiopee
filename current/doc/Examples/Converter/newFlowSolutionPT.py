# - newFlowSolution (pyTree) -
import Converter.Internal as Internal

# Create a flow solution node
n = Internal.newFlowSolution(name='FlowSolution', gridLocation='Vertex'); Internal.printTree(n)
#>> ['FlowSolution',None,[1 son],'FlowSolution_t']
#>>    |_['GridLocation',array('Vertex',dtype='|S1'),[0 son],'GridLocation_t']

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newFlowSolution(name='FlowSolution', gridLocation='Vertex', parent=z)
