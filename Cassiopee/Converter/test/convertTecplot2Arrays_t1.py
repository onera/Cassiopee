# - convertFile2Arrays (array) -
import Converter
import KCore.test as test

A=Converter.convertFile2Arrays('tec1-small.dat','fmt_tp')
#test.testA(A,10) # too small for diffArrays?
B=Converter.convertFile2Arrays('tec1-small.plt','bin_tp')
#test.testA(B,11)
#test.testA(Converter.diffArrays(A,B),12)

A=Converter.convertFile2Arrays('tec2-grid.dat','fmt_tp')
#test.testA(A,20)
B=Converter.convertFile2Arrays('tec2-grid.plt','bin_tp')
#test.testA(B,21)
#test.testA(Converter.diffArrays(A,B),22)

A=Converter.convertFile2Arrays('tec2-pressure.dat','fmt_tp')
test.testA(A,23)
B=Converter.convertFile2Arrays('tec2-pressure.plt','bin_tp')
test.testA(B,24)
#test.testA(Converter.diffArrays(A,B),25)

A=Converter.convertFile2Arrays('tec3-i.dat','fmt_tp')
test.testA(A,30)
B=Converter.convertFile2Arrays('tec3-i.plt','bin_tp')
test.testA(B,31)
#test.testA(Converter.diffArrays(A,B),32)

A=Converter.convertFile2Arrays('tec4-ij.dat','fmt_tp')
test.testA(A,40)
B=Converter.convertFile2Arrays('tec4-ij.plt','bin_tp')
test.testA(B,41)
#test.testA(Converter.diffArrays(A,B),42)

A=Converter.convertFile2Arrays('tec5-ijk.dat','fmt_tp')
test.testA(A,50)
B=Converter.convertFile2Arrays('tec5-ijk.plt','bin_tp')
test.testA(B,51)
#test.testA(Converter.diffArrays(A,B),52)

A=Converter.convertFile2Arrays('tec6-multiXY.dat','fmt_tp')
#test.testA(A,60)
#B=Converter.convertFile2Arrays('tec6-multiXY.plt','bin_tp')
#test.testA(B,61)
#test.testA(Converter.diffArrays(A,B),62)

# Not yet compatible
#A=Converter.convertFile2Arrays('tec7-multiXYShared.dat','fmt_tp')
#test.testA(A,70)
#B=Converter.convertFile2Arrays('tec7-multiXYShared.plt','bin_tp')
#test.testA(B,71)
#test.testA(Converter.diffArrays(A,B),72)

A=Converter.convertFile2Arrays('tec8-IJCellCenter.dat','fmt_tp')
test.testA(A,80)
B=Converter.convertFile2Arrays('tec8-IJCellCenter.plt','bin_tp')
test.testA(B,81)
#test.testA(Converter.diffArrays(A,B),82)

A=Converter.convertFile2Arrays('tec9-2D.dat','fmt_tp')
test.testA(A,90)
B=Converter.convertFile2Arrays('tec9-2D.plt','bin_tp')
test.testA(B,91)
#test.testA(Converter.diffArrays(A,B),92)

A=Converter.convertFile2Arrays('tec10-noheadersblock.dat','fmt_tp')
test.testA(A,100)
B=Converter.convertFile2Arrays('tec10-noheadersblock.plt','bin_tp')
test.testA(B,101)
#test.testA(Converter.diffArrays(A,B),102)

A=Converter.convertFile2Arrays('tec11-noheaderpoint.dat','fmt_tp')
test.testA(A,110)
B=Converter.convertFile2Arrays('tec11-noheaderpoint.plt','bin_tp')
test.testA(B,111)
#test.testA(Converter.diffArrays(A,B),112)

#A=Converter.convertFile2Arrays('tec12-noheaders.dat','fmt_tp')
#test.testA(A,120)
#B=Converter.convertFile2Arrays('tec12-noheaders.plt','bin_tp')
#test.testA(B,121)
#test.testA(Converter.diffArrays(A,B),122)

A=Converter.convertFile2Arrays('tec13-smallWithText.dat','fmt_tp')
test.testA(A,130)
B=Converter.convertFile2Arrays('tec13-smallWithText.plt','bin_tp')
test.testA(B,131)
#test.testA(Converter.diffArrays(A,B),132)
