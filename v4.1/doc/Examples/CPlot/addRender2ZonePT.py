# - addRender2Zone (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree as CPlot

a = G.cart((0,0,0), (1,1,1), (10,10,1))
C._initVars(a, '{Density}={CoordinateX}')
C._initVars(a, '{centers:VelocityX}={centers:CoordinateY}')
# With a material
a = CPlot.addRender2Zone(a, material='Glass', color='Blue', blending=1.,
                         meshOverlay=1, shaderParameters=[1.,1.])
# With field
a = CPlot.addRender2Zone(a, material='Solid', color='Iso:Density', blending=1.,
                         meshOverlay=1, shaderParameters=[1.,1.])
# With field+material
a = CPlot.addRender2Zone(a, material='Chrome', color='Iso:centers:VelocityX', blending=1.,
                         meshOverlay=1, shaderParameters=[1.,1.])

C.convertPyTree2File(a, 'out.cgns')
