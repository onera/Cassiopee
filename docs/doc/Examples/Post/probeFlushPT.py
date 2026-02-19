# - Probe flush (pyTree) -
import Post.Probe as Probe

# create a probe using mode 4
p4 = Probe.Probe('probe4.cgns', fields=['Density', 'Mach'], append=False, bufferSize=100)

for i in range(110):
    p4.extract(time=i*0.1, value=[i, i*10]) 
    # p4 will automatically flushed itself when reaching the buffer size limit (100)

# flush current buffered data to disk to save last extractions
p4.flush()