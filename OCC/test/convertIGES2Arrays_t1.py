# - convertIGES2Arrrays (IGES file format) -
import Converter as C
import OCC
import KCore.test as T

#~ coarse =  O.convertIGES2Arrays("hammer.iges", h=500., chordal_err=8.)
#~ C.convertArrays2File(coarse, 'hammer_coarsec.plt')

A = OCC.convertIGES2Arrays("hammer.iges", h=0., chordal_err=0.)
T.testA(A,1)

#~ fine = OCC.convertIGES2Arrays("hammer.iges", h=20., chordal_err=0.4)
#~ C.convertArrays2File(fine, 'hammer_fine.plt')
