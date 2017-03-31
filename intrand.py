import numpy as np
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

# q1 = cumtrapz(1*np.random.randn(10000),initial=0)
# q2 = cumtrapz(np.sqrt(2)*np.random.randn(10000),initial=0)
# q4 = cumtrapz(np.sqrt(4)*np.random.randn(10000),initial=0)
# q8 = cumtrapz(np.sqrt(8)*np.random.randn(10000),initial=0)

fig,(ax1,ax2,ax3) = plt.subplots(3,1)

ax1.plot(q1,label=r'x(t,a=1)')
ax1.plot(q2,label=r'x(t,a=2)')
ax1.plot(q4,label=r'x(t,a=4)')
ax1.plot(q8,label=r'x(t,a=8)')

ax1.set_title(r'$\frac{\partial x}{\partial t} = r,\quad r \sim \mathcal{N}(0,a)$')

ax1.legend(loc=3,ncol=4)

k = np.fft.fftfreq(10000)

N = 1000

# q1h = np.zeros(10000)
# q2h = np.zeros(10000)
# q4h = np.zeros(10000)
# q8h = np.zeros(10000)
# 
# for n in range(N):
#     
#     q1 = cumtrapz(1*np.random.randn(10000),initial=0)
#     q2 = cumtrapz(np.sqrt(2)*np.random.randn(10000),initial=0)
#     q4 = cumtrapz(np.sqrt(4)*np.random.randn(10000),initial=0)
#     q8 = cumtrapz(np.sqrt(8)*np.random.randn(10000),initial=0)
# 
#     q1h += abs(np.fft.fft(q1))**2
#     q2h += abs(np.fft.fft(q2))**2
#     q4h += abs(np.fft.fft(q4))**2
#     q8h += abs(np.fft.fft(q8))**2
# 
# q1h = q1h / N
# q2h = q2h / N
# q4h = q4h / N
# q8h = q8h / N

ax2.loglog(k[1:int(N/2-1)],q1h[1:int(N/2-1)],label=r'$|\hat{x}(t,a=1)|^2$')
ax2.loglog(k[1:int(N/2-1)],q2h[1:int(N/2-1)],label=r'$|\hat{x}(t,a=2)|^2$')
ax2.loglog(k[1:int(N/2-1)],q4h[1:int(N/2-1)],label=r'$|\hat{x}(t,a=4)|^2$')
ax2.loglog(k[1:int(N/2-1)],q8h[1:int(N/2-1)],label=r'$|\hat{x}(t,a=8)|^2$')

ax2.legend(loc=3,ncol=4)

ax3.semilogx(k[1:int(N/2-1)],q2h[1:int(N/2-1)]/q1h[1:int(N/2-1)],'C1',label=r'$\frac{|\hat{x}(t,a=2)|^2}{|\hat{x}(t,a=1)|^2}$')
ax3.semilogx(k[1:int(N/2-1)],q4h[1:int(N/2-1)]/q1h[1:int(N/2-1)],'C2',label=r'$\frac{|\hat{x}(t,a=4)|^2}{|\hat{x}(t,a=1)|^2}$')
ax3.semilogx(k[1:int(N/2-1)],q8h[1:int(N/2-1)]/q1h[1:int(N/2-1)],'C3',label=r'$\frac{|\hat{x}(t,a=8)|^2}{|\hat{x}(t,a=1)|^2}$')

ax3.set_ylim(-2,9)
ax3.legend(loc=3,ncol=3)

plt.show()