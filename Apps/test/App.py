# - App -
from Apps.App import App
myApp = App()
myApp.set(fileIn="toto.cgns")
myApp.get("fileIn")
myApp.requires(["fileIn", "fileOut", "bcType"])
myApp.check()
