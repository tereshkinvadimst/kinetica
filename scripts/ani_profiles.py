import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import os

# 1. Параметры
TH = 370
TC = 270
L  = 1e-4

# 2. Поиск и сортировка файлов (чтобы 100 шел после 90, а не после 10)
files = sorted(glob.glob("*.dat"), key=lambda x: os.path.getmtime(x) if "_" not in x else int(''.join(filter(str.isdigit, x))))

if not files:
    print("Файлы .dat не найдены!")
    exit()

# Подготовка окна
fig, ax = plt.subplots(figsize=(12, 6))
line_dsmc, = ax.plot([], [], label='DSMC расчёт', lw=2)
line_exact, = ax.plot([], [], '--', label='Аналитическое решение стационарного уравнения теплопроводности', alpha=0.8)

def init():
    ax.set_xlabel('X, [м]')
    ax.set_ylabel('T, [K]')
    ax.grid(True)
    ax.legend(loc='upper right')
    return line_dsmc, line_exact

def update(frame_file):
    # Загрузка данных текущего кадра
    data = np.genfromtxt(frame_file, delimiter='\t', skip_header=True, dtype=np.float64)
    
    x = data[:, 0]
    temp_dsmc = data[:, 11]
    temp_exact = TC + (TH - TC) / L * x
    
    # Обновление графиков
    line_dsmc.set_data(x, temp_dsmc)
    line_exact.set_data(x, temp_exact)
    
    # Динамическая подстройка осей
    ax.set_xlim(0, L)
    ax.set_ylim(260, 380)
    ax.set_title(f"Файл: {frame_file}")
    
    return line_dsmc, line_exact

# Создание анимации
# frames=files передает имя файла в функцию update
ani = FuncAnimation(fig, update, frames=files, init_func=init, blit=True, interval=200)

plt.tight_layout()
plt.show()

# # Если нужно сохранить:
# ani.save('kinetics_animation_Lem5.gif')
