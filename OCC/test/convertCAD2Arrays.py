# - convertCAD2Arrays (arrays) -
import Converter as C
import OCC

# IGES avec T3Mesher
A = OCC.convertCAD2Arrays('hammer.iges', format='fmt_iges',
                          h=0., chordal_err=0., growth_ratio=0.8, algo=1)
C.convertArrays2File(A, 'hammer1.plt')

# IGES avec OCC
A = OCC.convertCAD2Arrays('hammer.iges', format='fmt_iges',
                          chordal_err=1, algo=0)
C.convertArrays2File(A, 'hammer2.plt')

# STEP avec T3Mesher
A = OCC.convertCAD2Arrays('as1-oc-214.stp', format='fmt_step',
                          h=0., chordal_err=0., growth_ratio=0.8, algo=1)
C.convertArrays2File(A, 'as1.plt')

# STEP avec OCC
A = OCC.convertCAD2Arrays('as1-oc-214.stp', format='fmt_step',
                          chordal_err=1, algo=0)
C.convertArrays2File(A, 'as2.plt')
