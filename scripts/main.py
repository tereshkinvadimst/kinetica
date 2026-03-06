import kinetica as kt

MOLECULE_MASS = 6.63e-26 # Масса молекулы [кг]
FN            = 1e6      # Число реальных молекул на одну модельную частицу
MOLECULE_SIZE = 3.6e-10  # Размер атома аргона [м]
N_DENSITY     = 2.4e24   # Числовая плотность молекул [м^-3]
T_INI         = 273      # Начальная температура системы [K]
LX            = 5e-6     # Размер расчётной области по x [м]
LY            = 5e-6     # Размер расчётной области по y [м]
LZ            = 5e-6     # Размер расчётной области по z [м]
LCX           = 5e-9     # Характерный размер ячейки по x [м]
LCY           = 5e-7     # Характерный размер ячейки по y [м]
LCZ           = 5e-7     # Характерный размер ячейки по z [м]
NCX           = int(LX / LCX)                       # Число ячеек по x
NCY           = int(LY / LCY)                       # Число ячеек по y
NCZ           = int(LZ / LCZ)                       # Число ячеек по z
DT            = 1e-13    # Временной шаг моделирования [с]
TOTAL_TIME    = 5       # Полное время моделирования [с]
TC            = 263     # Температура холодной стенки
TH            = 283      # Температура горячей стенки

def main():
    domain_box = kt.Box(0, 0, 0, LX, LY, LZ)
    domain = kt.Domain(domain_box, MOLECULE_MASS, FN, MOLECULE_SIZE)
    domain.generateParticles(N_DENSITY, T_INI)
    #domain.saveXYZ("init.xyz")
    print(f"{NCX}, {NCY}, {NCZ}")
    domain.generateMesh(LCX, LCY, LCZ)
    domain.makeCellList()
    domain.setDiffuseWall(0, TH);
    domain.setDiffuseWall(1, TC);
    #domain.setDiffuseWall(2, T_INI);
    #domain.setDiffuseWall(3, T_INI);
    
    time: float = 0
    counter = 0
    while(time < TOTAL_TIME):
        domain.moveParticles(DT)
        domain.applyPeriodicBoundaries(True, True ,False)
        domain.updateCellList()
        domain.collideParticles(DT)
        if(counter % 100000 == 0):
            domain.printStatsHeader()
            domain.computeFlowProperties()
            domain.writeVTU(f"kinetica_{counter}.vtu")
        if(counter % 10000 == 0):
            domain.printStats(time)
           
        counter += 1
        time += DT
        
    
if __name__ == "__main__":
    main()