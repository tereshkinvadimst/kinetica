import kinetica as kt

MOLECULE_MASS = 6.63e-26 # Масса молекулы [кг]
FN            = 1e6      # Число реальных молекул на одну модельную частицу
MOLECULE_SIZE = 3.6e-10  # Размер атома аргона [м]
N_DENSITY     = 2.4e24   # Числовая плотность молекул [м^-3]
T_INI         = 300      # Начальная температура системы [K]
LX            = 1e-4     # Размер расчётной области по x [м]
LY            = 1e-4     # Размер расчётной области по y [м]
LZ            = 1e-6     # Размер расчётной области по z [м]
LCX           = 1e-5     # Характерный размер ячейки по x [м]
LCY           = 1e-5     # Характерный размер ячейки по y [м]
LCZ           = 1e-5     # Характерный размер ячейки по z [м]
NCX           = int(LX / LCX)                       # Число ячеек по x
NCY           = int(LY / LCY)                       # Число ячеек по y
NCZ           = int(LZ / LCZ)                       # Число ячеек по z
DT            = 1e-8     # Временной шаг моделирования [с]
TOTAL_TIME    = 1e-6     # Полное время моделирования [с]

def main():
    domain_box = kt.Box(0, 0, 0, LX, LY, LZ)
    domain = kt.Domain(domain_box, MOLECULE_MASS, FN, MOLECULE_SIZE)
    domain.generateParticles(N_DENSITY, T_INI)
    domain.saveXYZ("init.xyz")
    print(f"{NCX}, {NCY}, {NCZ}")
    domain.generateMesh(LCX, LCY, LCZ)
    domain.makeCellList()
    time: float = 0
    counter = 0
    while(time < TOTAL_TIME):
        domain.applyPeriodicBoundaries(True, True ,True)
        domain.updateCellList()
        domain.collideParticles(DT)
        domain.moveParticles(DT)
        domain.saveXYZ(f"kinetica_{counter}.xyz")
        print(f"time = {time}, counter = {counter}")
        counter += 1
        
    
if __name__ == "__main__":
    main()