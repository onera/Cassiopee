# - convertArrays2File (wav) -
import Converter as C
import math

# Sampling time step
Deltat = 5.e-5

# Number of samples
N = 500000

# Length of sound
L = N*Deltat
print('Total sound duration: %f'%L)

# Sound frequency (Hertz)
# Human ear is between 20 Hz and 20 kHz
f1 = 1000; f2 = 1100; f3 = 1200

# Signal
def F(time):
    return np.where(
        time < L/3.,
        np.cos(2*np.pi*f1*time),
        np.where(
            time < 2.*L/3.,
            np.cos(2*np.pi*f2*time),
            np.cos(2*np.pi*f3*time)
        )
    )

a = C.array('Time, Pressure', N, 1, 1)

# Time
for i in range(N): a[1][0,i] = Deltat*i

# Pressure
a = C.initVars(a, 'Pressure', F, ['Time'], isVectorized=True)

# Convert in wav uses Time and Pressure field
C.convertArrays2File([a], 'out.wav', 'bin_wav')
