module kinetica_particles
    use kinetica_kinds
    implicit none
    private

    public::generate_particles
    public::kill_particle
    
    !> Число живых частиц
    integer(i32),protected,public::n_particles = 0
    !> Запас частиц
    integer(i32),protected,public::capacity    = 0
    !> Положение частиц (3, Npart)
    real(f64),allocatable,public::r(:,:)
    !> Скорость частиц (3, Npart)
    real(f64),allocatable,public::u(:,:)
    !> Тип частиц (Npart)
    integer(i32),allocatable,public::species(:)
    !> Указатель на ячейку, которой принадлежит частица
    integer(i32),allocatable,public::cell_id(:)
    !> Указатель на следующую частицу в списке ячеек
    integer(i32),allocatable,public::next(:)
    !> Указатель на предыдущую частицу в списке ячеек
    integer(i32),allocatable,public::prev(:)

contains

    !> Сгенерировать частицы
    subroutine generate_particles(weight, n_density, temperature, u_flow, s, pos_min, pos_max, dt)
        use kinetica_random, only: poisson, gauss
        use kinetica_constants, only: r_gas
        use kinetica_species, only: molar_mass
        !> Вес частиц
        real(f64),intent(in),value::weight
        !> Числовая плотность частиц (м-3)
        real(f64),intent(in),value::n_density
        !> Температура
        real(f64),intent(in),value::temperature
        !> Скорость потока
        real(f64),intent(in)::u_flow(3)
        !> Тип частиц, которые нужно сгенерировать
        integer(i32),intent(in),value::s
        !> Минимальное положение области генерации частиц
        real(f64),intent(in)::pos_min(3)
        !> Максимально положение области генерации частиц
        real(f64),intent(in)::pos_max(3)
        !> Шаг по времени (опционален в случае перемещения частиц)
        real(f64),intent(in),optional,value::dt
        !> размеры области
        real(f64)::l(3)
        !> Объём области
        real(f64)::volume
        !> Приблизительное число частиц 
        real(f64)::n_nearly
        !> Случайное число на отрезке [0, 1)
        real(f64)::rnumber
        !> Число частиц, которые нужно сгенерировать
        integer(i32)::n_generated
        !> Среднеквадратичная скорость
        real(f64)::u_avg
        !> номер частицы
        integer(i32)::p
        !> случайно сгенерированное положение частицы
        real(f64)::position(3)
        !> случайно сгенерированная скорость частицы
        real(f64)::velocity(3)
        ! Вычисляем размер области
        l           = abs(pos_max - pos_min)
        ! Вычисляем объём области
        volume      = product(l)
        ! Вычисляем приблизительное число частиц в области
        n_nearly    = n_density * volume / weight
        if(n_nearly > 20) then ! Если частиц много, то:
            n_generated = int(n_nearly, i32)
            call random_number(rnumber)
            if(rnumber < n_nearly - real(n_generated, f64)) n_generated = n_generated + 1
        else ! Вычисляем число частиц из распределения Пуассона:
        n_generated = poisson(n_nearly)
        endif
        ! Вычисляем среднеквадратичную скорость
        u_avg = sqrt(3._f64 * r_gas * temperature / molar_mass(s))
    
        do p = 1, n_generated
            ! Генерируем равнораспределённое положение частицы
            call random_number(position)
            position = position * l + pos_min
            ! Генерируем скорость из распределения Гаусса
            velocity = u_flow + u_avg * gauss()
            ! Если указан шаг по времени, то также перемещаем частицу
            if(present(dt)) position = position + velocity * dt
            ! Добавляем частицу в конец
            call push_back_particle(position, velocity, s)
        end do

    end subroutine generate_particles


    subroutine push_back_particle(position, velocity, s)
        use kinetica_utils, only: resize
        use kinetica_species, only: n_particles_species
        real(f64),intent(in)::position(3)
        real(f64),intent(in)::velocity(3)
        integer(i32),intent(in),value::s
        integer(i32)::new_size
        integer(i32)::new_capacity

        new_size     = n_particles + 1
        new_capacity = capacity

        if(new_size > capacity) then
            ! Вычисляем новую ёмкость
            new_capacity = nint(new_size * 1.5_f64, kind=i32)
            ! Увеличиваем размер массивов
            call resize(r, [3, new_capacity])
            call resize(u, [3, new_capacity])
            call resize(species, new_capacity)
            call resize(cell_id, new_capacity)
            call resize(next, new_capacity)
            call resize(prev, new_capacity)

        endif

        r(:, new_size)         = position
        u(:, new_size)         = velocity
        species(new_size)      = s
        next(new_size)         = 0
        prev(new_size)         = 0
        n_particles            = new_size
        capacity               = new_capacity
        n_particles_species(s) = n_particles_species(s) + 1
    end subroutine push_back_particle

    subroutine kill_particle(i)
        use kinetica_species, only: n_particles_species
        !> Указатель на частицу, которую нужно удалить
        integer(i32),intent(in),value::i
        !> Указатель на последнюю частицу в списке частиц
        integer(i32)::last
        !> Указатели на соседей в двусвязном списке
        integer(i32)::p_i, n_i, p_last, n_last

        ! Если указатель невалидный, то пропускаем
        if (i < 1 .or. i > n_particles) return

        last = n_particles
        ! Уменьшаем число частиц на 1
        n_particles                     = n_particles - 1
        n_particles_species(species(i)) = n_particles_species(species(i)) - 1

        p_i    = prev(i)
        n_i    = next(i)
        p_last = prev(last)
        n_last = next(last)

        ! Удаляем узел i из списка
        if (p_i /= 0) next(p_i) = n_i
        if (n_i /= 0) prev(n_i) = p_i

        ! Если указатель не на последний элемент списка частиц
        if(i /= last) then
            ! Делаем swap
            r(:, i)    = r(:, last)
            u(:, i)    = u(:, last)
            species(i) = species(last)
            cell_id(i) = cell_id(last)
            ! Узел i должен принять связи старого last
            prev(i) = p_last
            next(i) = n_last

            ! Переподключаем соседей, которые раньше ссылались на last
            if (p_last /= 0) next(p_last) = i
            if (n_last /= 0) prev(n_last) = i
        endif
    end subroutine kill_particle
    
end module kinetica_particles