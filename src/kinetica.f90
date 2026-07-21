program kinetica
    use kinetica_kinds
    use kinetica_random
    use kinetica_particles
    use kinetica_species
    use kinetica_xyz
    use kinetica_domain
    use kineitca_mover
    use kinetica_piston
    use kinetica_init
    use kinetica_global
    implicit none
    real(f64),parameter::initial_density     = 1e27_f64
    real(f64),parameter::initial_temperature = 300._f64
    real(f64),parameter::initial_velocity(3) = [0._f64, 0._f64, 0._f64]
    real(f64),parameter::boundary_velocity(3) = [160._f64, 0._f64, 0._f64]
    integer(i32),parameter::n_steps          = 100
    integer(i32)::step

    call random_init(repeatable=.false., image_distinct=.true.)

    box_min_position = [0._f64, 0._f64, 0._f64]
    box_max_position = [1000e-9_f64, 100e-9_f64, 100e-9_f64]
    box_size         = box_max_position - box_min_position
    periodic         = [.false., .false., .true.]
    is_wall          = .true.
    is_wall(i_left)  = .false.
    is_wall(i_right) = .false.

    dt_coeff = 0.3_f64
    dl_coeff = 0.3_f64

    ! call set_piston(1, 0._f64, 160._f64, 1_i32)

    ! Создаём компонент
    call add_species("Ar", 0.03995_f64, 0.81_f64, 1.4_f64, 4.17e-10_f64, 273._f64)
    ! Генерируем частица аргона по всей области
    call generate_particles(1000._f64, initial_density, initial_temperature, initial_velocity, 1, box_min_position, box_max_position)
    ! Инициализируем параметры расчёта
    call compute_initial_params(initial_density, initial_temperature, max(maxval(initial_velocity), maxval(boundary_velocity), piston_velocity))
    ! Строим начальную сетку для вычисляения коллизий
    ! Расчётный цикл
    do step = 1, n_steps
        ! Генерируем частицы на входе
        call generate_particles(1000._f64, initial_density, initial_temperature, boundary_velocity, 1, [box_min_position(1) - dl, box_min_position(2), box_min_position(3)], [dl, box_max_position(2), box_max_position(3)])
        ! Двигаем частицы
        call move_particles(dt)
        ! Вычисляем коллизии
        ! Вычисляем новый шаг по времени и размер ячеек

        ! Двигаем поршень
        ! if(piston_on) call move_piston(dt)


        call write_xyz(to_string_i(step))
    end do

contains

pure function to_string_i(i) result(str)
    integer(i32), intent(in) :: i
    character(len=:), allocatable :: str

    character(len=32) :: tmp

    write(tmp,'(I0)') i
    str = trim(tmp)

end function

end program kinetica
