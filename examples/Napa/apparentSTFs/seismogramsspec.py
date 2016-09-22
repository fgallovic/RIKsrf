from numpy import loadtxt,savetxt,arange,fft,sqrt,vstack
from math import pi

seismograms=loadtxt('seismograms.dat')

dt=seismograms[1,0]
N=seismograms[:,0].size
T=seismograms[-1,0]
df=1./T
f=arange(N/2-1)*df
sf=abs(fft.fft(seismograms[:,1:],axis=0))*dt
sfspecD=sf[:N/2-1,:].T
sfspecA=(2*pi*f)**2*sfspecD
savetxt('seismogramsspecD.dat',vstack([f,sfspecD]).T)
savetxt('seismogramsspecA.dat',vstack([f,sfspecA]).T)

