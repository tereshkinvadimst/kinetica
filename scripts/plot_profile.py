import numpy as np
import matplotlib.pyplot as plt

fname = "kinetica_30000.dat"

data = np.genfromtxt(fname, delimiter='\t', skip_header=True, dtype=np.float64)

TH = 370
TC = 270
L  = 1e-4
TEXACT = TC + (TH - TC) / L * data[:, 0]

print(f"ux_mean = {np.mean(data[:, 2])}")
print(f"uy_mean = {np.mean(data[:, 3])}")
print(f"uz_mean = {np.mean(data[:, 4])}")

plt.figure(figsize=(14,7))
plt.plot(data[:, 0], data[:, 11], label='DSMC расчёт')
plt.plot(data[:, 0], TEXACT, label="Решение стационарного ур. теплопроводности")

plt.legend()
plt.xlabel('X, [м]')
plt.ylabel('T, [K]')
plt.grid(True)
plt.tight_layout()
plt.show()
