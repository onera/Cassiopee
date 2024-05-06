# - run tkCassiopee from a python script -
import CPlot.Tk as CTK
import tkCassiopee as K

t = C.convertFile2PyTree('out.cgns')

# Recuperation d'un arbre en memoire
CTK.t = t
(CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
CTK.display(CTK.t)

(win, frames, menu, menus, file, tools) = CTK.minimal2('Cassiopee ')
K.auto = {}
for app in K.ALLAPPS: K.auto[app] = 0
for app in K.TREEAPPS: K.addMenuItem(app, menus[0], frames[0])
for app in K.STATEAPPS: K.addMenuItem(app, menus[1], frames[1])
for app in K.EDGEAPPS: K.addMenuItem(app, menus[2], frames[2])
for app in K.SURFAPPS: K.addMenuItem(app, menus[3], frames[3])
for app in K.MESHAPPS: K.addMenuItem(app, menus[4], frames[4])
for app in K.BLOCKAPPS: K.addMenuItem(app, menus[5], frames[5])
for app in K.BCAPPS: K.addMenuItem(app, menus[6], frames[6])
for app in K.MOTIONAPPS: K.addMenuItem(app, menus[7], frames[7])
for app in K.SOLVERAPPS: K.addMenuItem(app, menus[8], frames[8])
for app in K.POSTAPPS: K.addMenuItem(app, menus[9], frames[9])
for app in K.VISUAPPS: K.addMenuItem(app, menus[10], frames[10])
for app in K.RENDERAPPS: K.addMenuItem(app, menus[11], frames[11])
    
win.mainloop()
CTK.PHOTOS = []
