import matplotlib.pyplot as plt
import numpy as np
import math

# Inisiasi nilai konstanta dan variabel dalam perhitungan
kb = 1      # Nilainya diambil 1 untuk memudahkan
N = 64      # Jumlah monomer dari urutan asam amino protein
eHH = 2     # Nilainya diambil 1 untuk memudahkan

datados = np.loadtxt('dos-reversed.txt')       # Me-load data file dos dari hasil simulasi wang-landau
datarg = np.loadtxt('Rg-reversed.txt')       # Me-load data file dos dari hasil simulasi wang-landau
datatau = np.loadtxt('Tau-reversed.txt')       # Me-load data file dos dari hasil simulasi wang-landau
E = -datados[:,0]      # Mengambil kolom pertama dari data file dos menjadi E (dikalikan negatif)
u = datados[:,1]       # Mengambil kolom kedua dari data file dos menjadi nilai dos
r = datarg[:,1]       # Mengambil kolom kedua dari data file dos menjadi nilai dos
tau = datatau[:,1]       # Mengambil kolom kedua dari data file dos menjadi nilai dos
a = 0.05
b = 2.4
c = 60
T = np.linspace(a,b,c)     # Menghasilkan nilai temperatur yang ingin diambil
Z = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Z dengan ukuran sama dengan T
Erata = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Erata dengan ukuran sama dengan T
Rgrata = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Erata dengan ukuran sama dengan T
Taurata = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Erata dengan ukuran sama dengan T
Ekuadrat_rata = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Ekuadrat_rata dengan ukuran sama dengan T
Cv = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Cv dengan ukuran sama dengan T
# print(T)

# A. Menghitung Fungsi Partisi Z(T)

Ztot = 0
for j in range(len(T)):
    # Menghitung Z untuk satu nilai T tertentu
    for i in range(len(E)):
        Zsem = u[i]*(math.exp(-E[i]/(kb*T[j])))
        Ztot = Ztot+Zsem
    Z[j] = Ztot
    Ztot = 0
# print("Z= ", Z)


# B. Menghitung rerata energi <E> dan <E^2>
# B.1. Menghitung <E>
# Menghitung <E> untuk semua T yang diambil
Erata_upper_tot = 0
for j in range(len(T)):
    # Menghitung <E> untuk satu nilai T tertentu
    for i in range(len(E)):
        Erata_upper_sem = E[i]*u[i]*(math.exp(-E[i]/(kb*T[j])))
        Erata_upper_tot = Erata_upper_tot+Erata_upper_sem
    Erata[j] = Erata_upper_tot/Z[j]
    Erata_upper_tot = 0
print("Erata= ", Erata)

# B.2 Menghitung <E^2>
# Menghitung <E^2> untuk semua T yang diambil
Ekuadrat_rata_upper_tot = 0
for j in range(len(T)):
    # Menghitung <E^2> untuk satu nilai T tertentu
    for i in range(len(E)):
        Ekuadrat_rata_upper_sem = (E[i]**2)*u[i]*(math.exp(-E[i]/(kb*T[j])))
        Ekuadrat_rata_upper_tot = Ekuadrat_rata_upper_tot+Ekuadrat_rata_upper_sem
    Ekuadrat_rata[j] = Ekuadrat_rata_upper_tot/Z[j]
    Ekuadrat_rata_upper_tot = 0
# print("Ekuadrat_rata= ", Ekuadrat_rata)

# C. Menghitung rata-rata Rg dan Tau
# Menghitung rata-rata Rg
Rgrata_upper_tot = 0
for j in range(len(T)):
    # Menghitung <Rg> untuk satu nilai T tertentu
    for i in range(len(E)):
        Rgrata_upper_sem = r[i]*u[i]*(math.exp(-E[i]/(kb*T[j])))
        Rgrata_upper_tot = Rgrata_upper_tot+Rgrata_upper_sem
    Rgrata[j] = Rgrata_upper_tot/Z[j]
    Rgrata_upper_tot = 0
print("Rgrata= ", Rgrata)
# Menghitung rata-rata Tau
Taurata_upper_tot = 0
for j in range(len(T)):
    # Menghitung <Tau> untuk satu nilai T tertentu
    for i in range(len(E)):
        Taurata_upper_sem = tau[i]*u[i]*(math.exp(-E[i]/(kb*T[j])))
        Taurata_upper_tot = Taurata_upper_tot+Taurata_upper_sem
    Taurata[j] = Taurata_upper_tot/Z[j]
    Taurata_upper_tot = 0
print("Taurata= ", Taurata)


# C. Menghitung Cv
# Menghitung Cv untuk semua nilai T yang diambil
for j in range(len(T)):
    # Menghitung Cv untuk satu nilai T tertentu
    Cv[j] = (Ekuadrat_rata[j]-(Erata[j]**2))/(kb*(T[j]**2))
print("Cv= ", Cv)

print("T= ",T)
# D. Plot grafik Cv vs T


"""
Tx = [0.0330, 0.2548, 0.4310, 1.0492]
Zx = [0.0330, 0.2548, 0.4310, 1.0492]
Eratax = [0.0330, 0.2548, 0.4310, 1.0492]
Ekuadrat_ratax = [0.0330, 0.2548, 0.4310, 1.0492]
Cvx = [0.0330, 0.2548, 0.4310, 1.0492]
CvxN = [0.0330, 0.2548, 0.4310, 1.0492]

Ztot = 0
for j in range(len(Tx)):
    # Menghitung Z untuk satu nilai T tertentu
    for i in range(len(E)):
        Zsem = u[i]*(math.exp(-E[i]/(kb*Tx[j])))
        Ztot = Ztot+Zsem
    Zx[j] = Ztot
    Ztot = 0

Erata_upper_tot = 0
for j in range(len(Tx)):
    # Menghitung <E> untuk satu nilai T tertentu
    for i in range(len(E)):
        Erata_upper_sem = E[i]*u[i]*(math.exp(-E[i]/(kb*Tx[j])))
        Erata_upper_tot = Erata_upper_tot+Erata_upper_sem
    Eratax[j] = Erata_upper_tot/Zx[j]
    Erata_upper_tot = 0

Ekuadrat_rata_upper_tot = 0
for j in range(len(Tx)):
    # Menghitung <E^2> untuk satu nilai T tertentu
    for i in range(len(E)):
        Ekuadrat_rata_upper_sem = (E[i]**2)*u[i]*(math.exp(-E[i]/(kb*Tx[j])))
        Ekuadrat_rata_upper_tot = Ekuadrat_rata_upper_tot+Ekuadrat_rata_upper_sem
    Ekuadrat_ratax[j] = Ekuadrat_rata_upper_tot/Zx[j]
    Ekuadrat_rata_upper_tot = 0

for j in range(len(Tx)):
    # Menghitung Cv untuk satu nilai T tertentu
    Cvx[j] = (Ekuadrat_ratax[j]-(Eratax[j]**2))/(kb*(Tx[j]**2))

for j in range (len(Tx)):
    CvxN[j] = Cvx[j]/N


print("Tx= ", Tx)
print("Eratax= ", Eratax)
#print("Cvx= ", Cvx)
#print("CvxN= ", CvxN)
"""

# plt.plot((T),(Cv/N), color='black', linestyle='solid', linewidth=2, marker='.', markerfacecolor='red', markersize=12)
# plt.plot((T),(Erata/N), color='blue', linestyle='solid', linewidth=2, marker='.', markerfacecolor='red', markersize=12)
# plt.plot((T),(Rgrata/N), color='blue', linestyle='solid', linewidth=2, marker='.', markerfacecolor='red', markersize=12)
# plt.plot((T),(Taurata), color='red', linestyle='solid', linewidth=2, marker='.', markerfacecolor='red', markersize=12)

fig, ax1 = plt.subplots()

ax1.plot(T/eHH, (Cv/N), color = "blue", linestyle='solid', linewidth=2, marker='.', markerfacecolor='blue', markersize=8)

ax2 = ax1.twinx()
ax2.plot(T/eHH, (Rgrata/N), color = "green", linestyle='solid', linewidth=2, marker='.', markerfacecolor='green', markersize=8)

ax3 = ax1.twinx()
ax3.plot(T/eHH, Taurata, color = "red", linestyle='solid', linewidth=2, marker='.', markerfacecolor='red', markersize=8)
ax3.spines["right"].set_position(("outward", 60))

ax1.set_ylabel("Cv/N", color = "blue")
ax2.set_ylabel("Rg/N", color = "green")
ax3.set_ylabel("Tortuosity", color = "red")

ax1.tick_params(axis="y", colors= "blue")
ax2.tick_params(axis="y", colors= "green")
ax3.tick_params(axis="y", colors= "red")

# Menambahkan beberapa detail pada plot
plt.xlim(0.0,1.2)

plt.show()

fig.savefig("plot3gambar.png", bbox_inches="tight")