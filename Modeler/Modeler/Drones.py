"""All drone models."""
import Geom as D
import Transform as T
import Generator as G
import Converter as C

#==============================================================================
# Bras
# IN: section: section dans le plan x,y
# IN: L: longeur du bras
# #==============================================================================
def droneArm(section, L):
    P0 = (0,0,0)
    P1 = (0,0,L)
    line = D.line(P0, P1)
    p = D.lineDrive(section, line)
    return p
