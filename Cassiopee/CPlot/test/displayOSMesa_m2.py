# - display (array) -
# display offscreen 360 using OSMesa
import CPlot.PyTree as CPlot
import Geom.PyTree as D
import Converter.Mpi as Cmpi
import KCore.test as test

LOCAL = test.getLocal()

if Cmpi.rank == 0: a = D.sphere((0,0,0), 1, N=200)
else: a = D.sphere((0,1,0), 0.8, N=200)

# 360 export
CPlot.display360(a, mode=0, stereo=0, offscreen=7, exportResolution='1024x512', export=LOCAL+'/out.png')
CPlot.finalizeExport(7)
if Cmpi.rank == 0: test.testF(LOCAL+'/out.png', 1)

# 360 WS
CPlot.display360(a, mode=0, stereo=2, offscreen=7, exportResolution='1024x512', export=LOCAL+'/out3.png')
CPlot.finalizeExport(7)
if Cmpi.rank == 0: test.testF(LOCAL+'/out3.png', 3)

# 360 ODS - trop lent
#CPlot.display360(a, mode=0, stereo=1, offscreen=7, exportResolution='256x128', export=LOCAL+'/out2.png')
#CPlot.finalizeExport(7)
#if Cmpi.rank == 0: test.testF(LOCAL+'/out2.png', 2)
