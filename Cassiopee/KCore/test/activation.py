import KCore
d = KCore.kcore.activation('0')
if d == 0: print('Open source version.')
else:
    year = d/12; month = d - 12*year
    print('Key is valid until %2d/%4d.'%(month, year))
