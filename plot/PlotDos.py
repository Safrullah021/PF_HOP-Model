import matplotlib.pyplot as plt
import numpy as np


data1 = np.loadtxt('dos-reversed.txt')
Emax1 = data1[:,1].max()

plt.scatter((-1*(data1[:,0])), np.log(data1[:,1]/Emax1),label='Simulasi Wang-Landau', color='red', marker='o', s=80)

plt.xlabel('Energi')
plt.ylabel('Log(g(E)/g(Emax))')
plt.legend()

plt.savefig("plotdos.png")

plt.show()
