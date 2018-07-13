# - ExtractDoublyDefined -
import Apps.Chimera.ExtractDoublyDefined as App

myApp = App.ExtractDoublyDefined()
myApp.set(dataIn='in.cgns', dataOut='visu.cgns', bcType='BCWall')
myApp.run()
