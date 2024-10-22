import matplotlib.pyplot as plt
import numpy as np

dataL = np.loadtxt('hzL.txt')
dataFr = np.loadtxt('hz.txt')

plt.plot(dataL[:,0], dataL[:,1], "r", label='$\Lambda CDM$')
plt.plot(dataFr[:,0], dataFr[:,1], 'b', label='$f(R)$')
plt.legend()
plt.ylabel('h(z) (km/s/Mpc)')
plt.xlabel('z')
plt.grid(True)
plt.savefig('hz.pdf')
plt.show()

dataL = np.loadtxt('dtdzL.txt')
dataFR = np.loadtxt('dtdz.txt')

plt.plot(dataL[:,1], dataL[:,0], 'r', label='$\Lambda CDM$')
plt.plot(dataFR[:,1], dataFR[:,0], 'b', label='$f(R)$')
plt.legend()
plt.ylabel('dtdz')
plt.xlabel('z')
plt.grid(True)
plt.savefig('dtdz.pdf')
plt.show()
