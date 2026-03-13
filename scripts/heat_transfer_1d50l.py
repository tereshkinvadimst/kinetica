import kinetica as kt
from math import pi, sqrt
import numpy as np

MOLECULE_MASS = 6.63e-26 # Масса молекулы [кг]
FN            = 1e2       # Число реальных молекул на одну модельную частицу
MOLECULE_SIZE = 3.6e-10  # Размер атома аргона [м]
N_DENSITY     = 2.4e24   # Числовая плотность молекул [м^-3]
T_INI         = 300      # Начальная температура системы [K]
SIGMA         = pi * MOLECULE_SIZE**2
LAMBDA        = 1 / sqrt(2) / N_DENSITY / SIGMA
LX            = 50 * LAMBDA # Размер расчётной области по x [м]
LY            = LAMBDA        # Размер расчётной области по y [м]
LZ            = LAMBDA        # Размер расчётной области по z [м]
LCX           = 5e-9        # Характерный размер ячейки по x [м]
LCY           = 5e-7        # Характерный размер ячейки по y [м]
LCZ           = 5e-7        # Характерный размер ячейки по z [м]
NCX           = int(LX / LCX)                       # Число ячеек по x
NCY           = int(LY / LCY)                       # Число ячеек по y
NCZ           = int(LZ / LCZ)                       # Число ячеек по z
DT            = 1e-13    # Временной шаг моделирования [с]
TOTAL_TIME    = 1e-7       # Полное время моделирования [с]
TC            = 200     # Температура холодной стенки
TH            = 500      # Температура горячей стенки
PISTON_VEL    = 200


def main():
    print(LAMBDA)
    domain_box = kt.Box(0., 0., 0., LX, LY, LZ)
    domain = kt.Domain(domain_box, MOLECULE_MASS, FN, MOLECULE_SIZE, 0.5)
    domain.generateParticles(N_DENSITY, T_INI, 0, LX)
    
    print(f"time step is {domain.getTimeStep()}")
    #domain.saveXYZ("init.xyz")
    domain.generateMesh()
    domain.makeCellList()
    
    wall1 = kt.DiffuseWall(
    np.array([LX, 0, 0], dtype=np.float64),np.array([-1, 0, 0], dtype=np.float64)
    , np.array([LY, LZ], dtype=np.float64), np.array([0, 0, 0], dtype=np.float64), TC
    )
    
    wall2 = kt.Wall(
    np.array([0, 0, 0], dtype=np.float64),np.array([1, 0, 0], dtype=np.float64)
    , np.array([LY, LZ], dtype=np.float64), np.array([PISTON_VEL, 0, 0], dtype=np.float64)
    )
    
    domain.addWall(wall1)
    domain.addWall(wall2)

    
    time: float  = 0
    counter: int = 0
    while(time < TOTAL_TIME):

        domain.moveParticles()
        domain.applyPeriodicBoundaries(False, True , True)
        domain.updateCellList()
        domain.collideParticles()

        if(counter % 100 == 0):
            domain.printStatsHeader()
            domain.computeFlowProperties()
            domain.writeXProfile(f"kinetica_{counter}.dat")
        if(counter % 10 == 0):
            domain.printStats(time)

           
        counter += 1
        time += domain.getTimeStep()
        
    
if __name__ == "__main__":
    main()