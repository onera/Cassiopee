# - Probe printInfo (pyTree) -
import Post.Probe as Probe

# create a probe using mode 4
p4 = Probe.Probe('probe4.cgns', fields=['Density', 'Mach'], append=False, bufferSize=100)

# print all relevant information about the p4 probe object
p4.printInfo()