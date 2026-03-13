import kinetica as kt

MOLECULE_MASS = 6.63e-26 # Масса молекулы [кг]
FN            = 1e2       # Число реальных молекул на одну модельную частицу
MOLECULE_SIZE = 3.6e-10  # Размер атома аргона [м]
N_DENSITY     = 2.4e24   # Числовая плотность молекул [м^-3]
T_INI         = 300      # Начальная температура системы [K]
LX            = 1e-4     # Размер расчётной области по x [м]
LY            = 1e-6     # Размер расчётной области по y [м]
LZ            = 1e-6     # Размер расчётной области по z [м]
LCX           = 5e-9     # Характерный размер ячейки по x [м]
LCY           = 5e-7     # Характерный размер ячейки по y [м]
LCZ           = 5e-7     # Характерный размер ячейки по z [м]
NCX           = int(LX / LCX)                       # Число ячеек по x
NCY           = int(LY / LCY)                       # Число ячеек по y
NCZ           = int(LZ / LCZ)                       # Число ячеек по z
DT            = 1e-13    # Временной шаг моделирования [с]
TOTAL_TIME    = 1       # Полное время моделирования [с]
TC            = 270     # Температура холодной стенки
TH            = 370      # Температура горячей стенки

N_DENSITY_L   = 2.4e24
N_DENSITY_R   = 2.4e24
T_L           = 4800
T_R           = 48

def main():
    domain_box = kt.Box(0., 0., 0., LX, LY, LZ)
    domain = kt.Domain(domain_box, MOLECULE_MASS, FN, MOLECULE_SIZE, 0.5)
    domain.generateParticles(N_DENSITY_L, T_L, 0, LX / 2)
    domain.generateParticles(N_DENSITY_R, T_R, LX / 2, LX)
    
    print(f"time step is {domain.getTimeStep()}")
    #domain.saveXYZ("init.xyz")
    domain.generateMesh()
    domain.makeCellList()

    domain.setDiffuseWall(0, T_L)
    domain.setDiffuseWall(1, T_R)
    
    time: float  = 0
    counter: int = 0
    while(time < TOTAL_TIME):

        domain.moveParticles()
        domain.applyPeriodicBoundaries(False, True , True)
        domain.updateCellList()
        domain.collideParticles()


        if(counter % 100 == 0):
            domain.printStatsHeader()
        if(counter % 10 == 0):
            domain.computeFlowProperties()
            domain.writeVTU(f"kinetica_{counter}.vtu")
            domain.printStats(time)
           
        counter += 1
        time += domain.getTimeStep()
        
    
if __name__ == "__main__":
    main()