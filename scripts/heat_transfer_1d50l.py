import kinetica as kt
from math import pi, sqrt
import numpy as np

MOLECULE_MASS = 6.63e-26 # Масса молекулы [кг]
FN            = 1e6       # Число реальных молекул на одну модельную частицу
MOLECULE_SIZE = 3.6e-10  # Размер атома аргона [м]
N_DENSITY     = 3.2e22   # Числовая плотность молекул [м^-3]
T_INI         = 300      # Начальная температура системы [K]
SIGMA         = pi * MOLECULE_SIZE**2
LAMBDA        = 1 / sqrt(2) / N_DENSITY / SIGMA
LX            = 200 * LAMBDA # Размер расчётной области по x [м]
LY            = LAMBDA        # Размер расчётной области по y [м]
LZ            = LAMBDA        # Размер расчётной области по z [м]
LCX           = 5e-9        # Характерный размер ячейки по x [м]
LCY           = 5e-7        # Характерный размер ячейки по y [м]
LCZ           = 5e-7        # Характерный размер ячейки по z [м]
NCX           = int(LX / LCX)                       # Число ячеек по x
NCY           = int(LY / LCY)                       # Число ячеек по y
NCZ           = int(LZ / LCZ)                       # Число ячеек по z
DT            = 1e-13    # Временной шаг моделирования [с]
TOTAL_TIME    = 5e-6       # Полное время моделирования [с]
TC            = 200     # Температура холодной стенки
TH            = 500      # Температура горячей стенки
PISTON_VEL    = 4481.36401988098


def main():
    print(LAMBDA)
    domain_box = kt.Box(0., 0., 0., LX, LY, LZ)
    domain = kt.Domain(domain_box, MOLECULE_MASS, FN, MOLECULE_SIZE, 0.4, 0.4)
    domain.generateParticles(N_DENSITY, T_INI, LAMBDA, LX)
    
    print(f"time step is {domain.getTimeStep()}")
    #domain.saveXYZ("init.xyz")
    domain.generateMesh()
    domain.makeCellList()
    
 
    # ===== Стенки =====
    walls = []

    # Поршень и противоположная стенка по X
    walls.append(kt.Wall(np.array([0, 0, 0], dtype=np.float64),
                                np.array([1, 0, 0], dtype=np.float64),
                                np.array([LY, LZ], dtype=np.float64),
                                np.array([PISTON_VEL, 0, 0], dtype=np.float64),
                                TH))  # поршень горячий
    walls.append(kt.Wall(np.array([LX, 0, 0], dtype=np.float64),
                                np.array([-1, 0, 0], dtype=np.float64),
                                np.array([LY, LZ], dtype=np.float64),
                                np.array([0, 0, 0], dtype=np.float64),
                                TC))  # противоположная стенка холодная

    # Стенки по Y
    walls.append(kt.Wall(np.array([0, 0, 0], dtype=np.float64),
                                np.array([0, 1, 0], dtype=np.float64),
                                np.array([LX, LZ], dtype=np.float64),
                                np.array([0, 0, 0], dtype=np.float64),
                                TC))
    walls.append(kt.Wall(np.array([0, LY, 0], dtype=np.float64),
                                np.array([0, -1, 0], dtype=np.float64),
                                np.array([LX, LZ], dtype=np.float64),
                                np.array([0, 0, 0], dtype=np.float64),
                                TC))

    # Стенки по Z
    walls.append(kt.Wall(np.array([0, 0, 0], dtype=np.float64),
                                np.array([0, 0, 1], dtype=np.float64),
                                np.array([LX, LY], dtype=np.float64),
                                np.array([0, 0, 0], dtype=np.float64),
                                TC))
    walls.append(kt.Wall(np.array([0, 0, LZ], dtype=np.float64),
                                np.array([0, 0, -1], dtype=np.float64),
                                np.array([LX, LY], dtype=np.float64),
                                np.array([0, 0, 0], dtype=np.float64),
                                TC))

    # Добавляем все стены в домен
    for w in walls:
        domain.addWall(w)
        
    time: float  = 0
    counter: int = 0
    while(time < TOTAL_TIME):
                    
        piston: kt.Wall = walls[0]  # первый в списке - поршень

        domain.moveParticles()
        domain.applyPeriodicBoundaries(False, False , False)
        domain.updateCellList()
        domain.collideParticles()

        if(counter % 10 == 0):
            domain.printStatsHeader()
            domain.computeFlowProperties()
            domain.writeXProfile(f"kinetica_{counter}.dat")
            print(f"Wall position: {piston.getCenterPosition()}")
        if(counter % 1 == 0):
            domain.printStats(time)
            
        # Остановка поршня
        if piston.getCenterPosition()[0] > LX / 4:
            piston.setVelocity(np.array([0, 0, 0]))
           
        counter += 1
        time += domain.getTimeStep()
        
    
if __name__ == "__main__":
    main()