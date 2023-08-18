import matplotlib.pyplot as plt
import numpy as np


data1 = np.loadtxt('dos-reversed.txt')
data2 = np.loadtxt('eksak-reversed.txt')
Emax1 = data1[:,1].max()
Emax2 = data2[:,1].max()

plt.scatter((-1*(data2[:,0])), np.log(data2[:,1]/Emax2),label='Hasil Eksak', color='black', marker='o', s=80)
plt.scatter((-1*(data1[:,0])), np.log(data1[:,1]/Emax1),label='Simulasi Wang-Landau', color='red', marker='x', s=80)

plt.xlabel('Energi')
plt.ylabel('Log(g(E)/g(Emax))')
plt.legend()

plt.savefig("plotdos.png")

plt.show()