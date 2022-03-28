# - CODA.IBM -
import Apps.Coda.ToolboxIBM_CODA as AppIBM
import Converter.PyTree as C
import Converter.Internal as Internal

#tb = C.convertFile2PyTree("case.cgns")
tb = C.convertFile2PyTree("naca1DRANS.cgns")
for z in Internal.getZones(tb):
    AppIBM._snearFactor(z,1.)

AppIBM.prepareOctree(tb, t_out='t.cgns',check=True, vmin=5)
