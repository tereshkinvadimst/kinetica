from math import sqrt
from scipy.constants import R

# Скорость ударной волны
U_SW = 6000
# Показатель адиабаты
GAMMA = 1.67
# Температура газа
T = 300
# Молярная масса газа
M = 0.039948
#-------------------------
# Скорость звука в газе
Cs = sqrt(GAMMA * R * T / M)
# Скорость поршня
Up = 2 / (GAMMA + 1) * (U_SW - Cs * Cs / U_SW)

print(f"Скорость поршня: {Up} м/c")