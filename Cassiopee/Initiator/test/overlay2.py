# - overlayField -
import Converter as C
import Initiator as I

from math import *
Mach = 0.626

list = ['naca.tp','cartL1.tp','cartL2.tp','cartL3.tp','cartL4.tp','cartL5.tp','cartL6.tp','cartL7.tp','cartL8.tp','cartL9.tp','cartL10.tp','cartL11.tp','cartL12.tp']
zone = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
count = 0

for x in list:

    count = count+1
    countName = '%(#)d'%{"#": count}
    at = C.convertFile2Arrays(x,"fmt_tp")
    a = at[0]

    ac = I.initScully(a, (-10.,-0.25), 0.2536, 0.162, Mach, 0)

    at = C.convertFile2Arrays("field.plt","bin_tp")
    a2 = at[count-1]

    an = I.overlayField(ac, a2, Mach)
    file = 'newL'+countName+'.plt'
    C.convertArrays2File([an],file,"bin_tp")
