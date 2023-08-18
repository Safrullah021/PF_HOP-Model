import matplotlib.pyplot as plt
import numpy as np
import math

# Inisiasi nilai konstanta dan variabel dalam perhitungan
kb = 1      # Nilainya diambil 1 untuk memudahkan
N = 62      # Jumlah monomer dari urutan asam amino protein
eHH = 2     # Nilainya diambil 1 untuk memudahkan

datados = np.loadtxt('dos-reversed.txt')       # Me-load data file dos dari hasil simulasi wang-landau


E = -datados[:,0]      # Mengambil kolom pertama dari data file dos menjadi E (dikalikan negatif)
u = datados[:,1]       # Mengambil kolom kedua dari data file dos menjadi nilai dos


a = 0.05
b = 2.4
c = 50
T = np.linspace(a,b,c)     # Menghasilkan nilai temperatur yang ingin diambil
Z = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Z dengan ukuran sama dengan T
Erata = np.linspace(a,b,c)     # Membuat variabel tipe matriks untuk Erata dengan ukuran sama dengan T
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

# C. Menghitung Cv
# Menghitung Cv untuk semua nilai T yang diambil
for j in range(len(T)):
    # Menghitung Cv untuk satu nilai T tertentu
    Cv[j] = (Ekuadrat_rata[j]-(Erata[j]**2))/(kb*(T[j]**2))
print("Cv= ", Cv)

print("T= ",T)
# D. Plot grafik Cv vs T

fig, ax1 = plt.subplots()

ax1.plot(T/eHH, (Cv/N), color = "blue", linestyle='solid', linewidth=2, marker='.', markerfacecolor='blue', markersize=8)


ax1.set_ylabel("Cv/N", color = "blue")
ax1.tick_params(axis="y", colors= "blue")

# Menambahkan beberapa detail pada plot
plt.xlim(0.0,1.2)

plt.show()

fig.savefig("plot3gambar.png", bbox_inches="tight")