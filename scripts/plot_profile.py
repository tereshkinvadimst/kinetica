import numpy as np
import matplotlib.pyplot as plt

fname = "kinetica_57000.dat"  # замените при необходимости

data = np.genfromtxt(fname, delimiter='\t', skip_header=True, dtype=np.float64)

print(np.mean(data[:, 1]))
plt.figure(figsize=(20,7))
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('X [m]')
plt.ylabel('Ttr [K]')
plt.title('Температура Ttr по координате X')
plt.grid(False)
plt.tight_layout()
plt.show()
