from numpy import loadtxt,savetxt,arange,fft,sqrt,vstack
from math import pi

stf=loadtxt('stf.dat')

dt=stf[1,0]
N=stf[:,0].size
T=stf[-1,0]
df=1./T
f=arange(N/2-1)*df
spec=abs(fft.fft(stf[:,1]))*dt
specD=spec[:N/2-1]
specA=(2*pi*f)**2*specD
savetxt('stfspecD.dat',vstack([f,specD]).T)
savetxt('stfspecA.dat',vstack([f,specA]).T)
