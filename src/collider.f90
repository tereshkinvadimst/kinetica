module kinetica_collider
    use kinetica_kinds
    implicit none
    private

    public::compute_collisions
    
contains

!> Функция вычисления сечения столкновения для пары частиц
pure function compute_sigmag(s1, s2, gij) result(sigmag)
    use kinetica_species, only: d => molecule_size             &
                              , omega                          &
                              , t_ref => reference_temperature &
                              , m => molar_mass
    use kinetica_constants, only: pi, n_a, k_b
    integer(i32),intent(in),value::s1 !< Тип первой частицы
    integer(i32),intent(in),value::s2 !< Тип второй частицы
    real(f64),intent(in),value::gij   !< Относительная скорость молекул
    real(f64)::sigmag                 !< Сечение столкновения молекул
    real(f64)::d_avg                  !< Средний опорный диаметр столкновения пар частиц
    real(f64)::omega_avg              !< Средний опорный показатель вязкости пары частиц
    real(f64)::t_ref_avg              !< Средняя опорная температура пары частиц
    real(f64)::mr                     !< Приведённая масса молекулы пары частиц
    real(f64)::hs_model               !< Часть уравнения, соответсвующая модели твёрдых шаров (HS)
    real(f64)::vhs_model              !< Часть уравнения, соответсвующая модели VHS

    ! Вычисляем параметры столкновения пар частиц
    d_avg     = 0.5_f64 * (d(s1) + d(s2))
    omega_avg = 0.5_f64 * (omega(s1) + omega(s2))
    t_ref_avg = 0.5_f64 * (t_ref(s1) + t_ref(s2))
    mr        = m(s1) * m(s2) / (m(s1) + m(s2)) / n_a

    ! Вклад модели HS
    hs_model = pi * d_avg**2
    ! Вклад модели VHS
    vhs_model = (2._f64 * k_b * t_ref_avg/mr)**(omega_avg - 0.5_f64) / gamma(2.5_f64 - omega_avg) &
              * gij**(2._f64 - 2._f64 * omega_avg)

     sigmag =  hs_model * vhs_model
end function compute_sigmag


!> Процедура вычисления относительной скорости
!> пары частиц после столкновения по модели VSS.
!> Отнормированная тройка строилась по методу:
!> Duff, T., Burgess, J., Christensen, P., Hery, C.,
!> Kensler, A., Liani, M., Villemin, R. 
!> «Building an Orthonormal Basis, Revisited». 
!> Journal of Computer Graphics Techniques (JCGT),
!> vol. 6, no. 1, 2017, pp. 52–59.
!> https://jcgt.org/published/0006/01/01/
subroutine compute_scatter(s1, s2, gij, g_in, g_out)
    use kinetica_species, only: alpha
    use kinetica_constants, only: pi
    integer(i32),intent(in),value::s1 !< Тип первой частицы
    integer(i32),intent(in),value::s2 !< Тип второй частицы
    real(f64),intent(in),value::gij   !< Модуль относительной скорости
    real(f64),intent(in)::g_in(3)     !< Вектор относительной скорости до столкновения
    real(f64),intent(out)::g_out(3)   !< Вектор относительной скорости после столкновения
    real(f64)::alpha_avg              !< Средний параметр рассеяния
    real(f64)::coschi                 !< Косинус угла отклонения
    real(f64)::sinchi                 !< Синус угла отклонения
    real(f64)::phi                    !< Случайная величина на отрезке [0, 2pi)
    real(f64)::r                      !< Случайная величина на отрезке [0, 1)
    real(f64)::s                      !< Знак компоненты
    real(f64)::a                      !< Промежуточная переменная
    real(f64)::b                      !< Промежуточная переменная
    real(f64)::e1(3)                  !< Нормаль полёта
    real(f64)::b1(3)                  !< Вектор, который задаёт плоскость поворота
    real(f64)::b2(3)                  !< Вектор, который задаёт плоскость поворота

    ! Вычисляем средний параметр рассеяния
    alpha_avg = 0.5_f64 * (alpha(s1) + alpha(s2))
    ! Генерируем случайное число на отрезке [0, 1)
    call random_number(r)
    ! Вычисляем косинус и синус угла рассеяния
    coschi = 2._f64 * r**(1._f64 / alpha_avg) - 1._f64
    sinchi = sqrt(max(0._f64, 1._f64 - coschi*coschi))
    ! Генерируем случайное число на отрезке [0, 2pi)
    call random_number(r)
    phi = 2.0_f64 * pi * r

    ! Вычисляем нормаль полёта
    e1 = g_in / gij

    ! Вычисляем два перпендикулярных вектора (https://jcgt.org/published/0006/01/01/)
    s  = sign(1._f64, e1(3))
    a  = -1._f64 / (s + e1(3))
    b  = e1(1) * e1(2) * a
    b1 = [1._f64 + s * e1(1)**2 * a, s * b           , -s * e1(1)]
    b2 = [b                        , s + e1(2)**2 * a, -e1(2)    ]

    ! Вычисляем поворот вектора после столкновения
    g_out = gij * (coschi * e1 + sinchi * (cos(phi) * b1 + sin(phi) * b2))
end subroutine compute_scatter

subroutine rand_pair(npc, i, j)
    integer(i32),intent(in)::npc !< Число частиц
    integer(i32),intent(out)::i  !< Случайный номер первой частицы
    integer(i32),intent(out)::j  !< Случайный номер второй частицы
    real(f64)::r                 !< Случайное число на отрезке [0, 1)
    ! Генерируем случайное число
    call random_number(r)
    ! Получаем случайный номер частицы
    i = min(1 + int(r * npc, i32), npc)
    ! Генерируем номер второй частицы
    do
        call random_number(r)
        j = min(1 + int(r * npc, i32), npc)
        if (j /= i) exit
    end do
end subroutine rand_pair

!> Вычислени коллизий между частицами в ячейке
subroutine compute_collisions(particles_idxs, volume, dt, weight, npc_avg, sigmag_max, n_coll)
    use kinetica_random, only: poisson
    use kinetica_particles, only: species, u
    use kinetica_species, only: m => molar_mass
    use kinetica_constants, only: n_a
    integer(i32),intent(in)::particles_idxs(:) !< Индексы частиц в ячейке
    real(f64),intent(in),value::volume         !< Объём ячейки
    real(f64),intent(in),value::dt             !< Шаг по времени
    real(f64),intent(in),value::weight         !< Вес модельных частиц
    real(f64),intent(in),value::npc_avg        !< Среднее число частиц в ячейке
    real(f64),intent(inout)::sigmag_max        !< Максимальное сечение столкновения
    integer(i32),intent(out)::n_coll           !< Число принятых столкновений
    integer(i32)::npc                          !< Число частиц в ячейке
    real(f64)::n_coll_avg                      !< Среднее число столкновений в ячейке
    integer(i32)::n_cand                       !< Число кандидатов для столкновения
    integer(i32)::ic                           !< Номер столкновения
    integer(i32)::i                            !< Номер первой частицы для столкновения
    integer(i32)::j                            !< Номер второй частицы для столкновения
    integer(i32)::k1                           !< Номер первой частицы в глобальном списке
    integer(i32)::k2                           !< Номер второй частицы в глобальном списке
    integer(i32)::s1                           !< Тип первой частицы
    integer(i32)::s2                           !< Тип второй частицы
    real(f64)::gijv(3)                         !< Вектор относительной скорости
    real(f64)::gijv_new(3)                     !< Вектор относительной скорости после столкновения
    real(f64)::gij                             !< Модуль вектора относительной скорости
    real(f64)::sigmag                          !< Сечение столкновения частиц
    real(f64)::r                               !< Случайное число на отрезке [0, 1)
    real(f64)::m1                              !< Масса первой молекулы
    real(f64)::m2                              !< Масса второй молекулы
    real(f64)::a                               !< Временная переменная для расчёта скорости центра масс
    real(f64)::b                               !< Временная переменная для расчёта скорости центра масс
    real(f64)::ucm(3)                          !< Скорость центра масс


    n_coll = 0

    npc = size(particles_idxs)

    ! Если в ячейке меньше двух частиц, то столкновения не считаем
    if(npc < 2) return

    ! Вычисляем среднее число столкновений
    n_coll_avg =  0.5_f64 * npc_avg * (npc - 1) * weight * sigmag_max * dt / volume
    ! Вычисляем число кандидатов из распределения Пуассона
    n_cand     = poisson(n_coll_avg)

    do ic = 1, n_cand
        ! Вычисляем случайную пару
        call rand_pair(npc, i, j)
        ! Получаем реальные индексы частиц
        k1 = particles_idxs(i)
        k2 = particles_idxs(j)
        ! Получаем типы частиц
        s1 = species(k1)
        s2 = species(k2)
        ! Вычисляем относительную скорость частиц
        gijv = u(:, k1) - u(:, k2)
        gij  = norm2(gijv)
        ! Вычисляем сечение столкновения
        sigmag = compute_sigmag(s1, s2, gij)
        ! Обновляем максимум
        sigmag_max = max(sigmag, sigmag_max)
        ! Генерируем случайное число на отрезке [0, 1)
        call random_number(r)
        ! Принимаем с вероятностью
        if(sigmag < sigmag_max * r) cycle
        ! Вычисляем относительную скорость после столкновения
        call compute_scatter(s1, s2, gij, gijv, gijv_new)
        ! Вычисляем массу молекул
        m1 = m(s1) / n_a
        m2 = m(s2) / n_a
        ! Вычисляем скорость центра масс
        a   = m1 / (m1 + m2)
        b   = m2 / (m1 + m2)
        ucm = a * u(:, k1) + b * u(:, k2)
        ! Вычисляем скорости после столкновения
        u(:, k1) = ucm + b * gijv_new
        u(:, k2) = ucm - a * gijv_new
        ! Увеличиваем счётчик принятых столкновений
        n_coll = n_coll + 1
    end do
end subroutine compute_collisions

end module kinetica_collider