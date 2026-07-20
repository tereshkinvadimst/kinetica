module kinetica_xyz
    use kinetica_kinds
    implicit none
    private

    public::write_xyz
    
contains

subroutine write_xyz(unique_name)
    use kinetica_particles, only: n_part => n_particles, r, u, species
    use kinetica_species,   only: species_names
    character(len=*),intent(in)::unique_name  !< базовое имя файла (без .xyz)
    integer::unit,ios                         !< дескриптор и код ошибки
    integer(i32)::i, s                        !< счётчик частиц и её тип
    character(len=:),allocatable::fname

    fname = trim(unique_name)//'.xyz'
    open(newunit=unit, file=fname, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        write(*,'(2a)') 'write_xyz: не удалось открыть файл ', fname
        return
    end if

    ! Заголовок кадра: число атомов + описание колонок для OVITO
    write(unit,'(i0)') n_part
    write(unit,'(a)')  'Properties=species:S:1:pos:R:3:velocity:R:3'

    ! Данные: символ элемента, координаты (3), скорость (3)
    do i = 1, n_part
        s = species(i)
        write(unit,'(a,6(1x,es16.8))') trim(species_names(s)), &
              r(1,i) / 1e-9_f64, r(2,i) / 1e-9_f64, r(3,i) / 1e-9_f64, u(1,i), u(2,i), u(3,i)
    end do

    close(unit)
end subroutine write_xyz
    
end module kinetica_xyz