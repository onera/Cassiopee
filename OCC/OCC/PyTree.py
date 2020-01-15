# OCC pyTree
try:
    import OCC as O
    import numpy
    import re
    import KCore
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
except: raise ImportError("OCC.PyTree: requires Converter module.")

__version__ = O.__version__

#==============================================================================
# -- convertCAD2PyTree --
#==============================================================================
def convertCAD2PyTree(fileName, format='fmt_iges', h=0., chordal_err=0., growth_ratio=0.,
                      deflection=1., algo=0):
  """Convert a CAD (IGES or STEP) file to pyTree.
  Usage: convertCAD2PyTree(fileName, options)"""

  a = O.convertCAD2Arrays(fileName, format, h, chordal_err, growth_ratio)
  
  t = C.newPyTree([])
  base1 = False; base2 = False; base3 = False; base = 1; c = 0
  
  for i in a:
    if len(i) == 5: # Structure
      if i[3] == 1 and i[4] == 1:
        if base1 == False:
          t = C.addBase2PyTree(t, 'Base1', 1); base1 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base1][2].append(z)
      elif i[4] == 1:
        if base2 == False:
          t = C.addBase2PyTree(t, 'Base2', 2); base2 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__) 
        t[2][base2][2].append(z)
      else:
        if base3 == False:
          t = C.addBase2PyTree(t, 'Base', 3); base3 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base3][2].append(z)
    else: # non structure
      if i[3] == 'BAR':
        if base1 == False:
          t = C.addBase2PyTree(t, 'Base1', 1); base1 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base1][2].append(z)
      elif i[3] == 'TRI' or i[3] == 'QUAD':
        if base2 == False:
          t = C.addBase2PyTree(t, 'Base2', 2); base2 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base2][2].append(z)
      else:
        if base3 == False:
          t = C.addBase2PyTree(t, 'Base', 3); base3 = base; base += 1
        z = Internal.createZoneNode(C.getZoneName('Zone'), i, [],
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
        t[2][base3][2].append(z)

  Internal._correctPyTree(t, level=2) # force unique name
  Internal._correctPyTree(t, level=7) # create familyNames
  return t
