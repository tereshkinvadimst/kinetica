module kinetica_init
    use kinetica_kinds
    implicit none
    private

    public::compute_initial_params

contains

!> Процедура вычисления начальных параметров для dsmc расчёта
subroutine compute_initial_params(n_density, temperature, velocity)
    use kinetica_species, only: n_species, molar_mass, molecule_size
    use kinetica_constants, only: r_gas, k_b, pi
    use kinetica_global
    real(f64),intent(in),value::n_density   !< Начальная числовая плотность молекул
    real(f64),intent(in),value::temperature !< Начальная температура
    real(f64),intent(in),value::velocity    !< Начальная скорость молекул/поршня
    real(f64)::min_mass                     !< Минимальная масса молекулы
    real(f64)::max_d                        !< Максимальный размер молекул
    real(f64)::rms_velocity                 !< Средняя квадратичная скорость молекул

    ! Вычисляем минимальную массу молекулы
    min_mass     = minval(molar_mass)
    ! Вычисляем максимальный размер молекулы
    max_d        = maxval(molecule_size)
    ! Вычисляем среднеквадратичную скорость молекул
    rms_velocity = sqrt(3._f64 * r_gas * temperature / min_mass)

    ! Максимальная скорость равняется среднекинетической + скорости потока 
    max_velocity = rms_velocity + velocity
    ! Вычисляем минимальную длину свободного пробега
    lambda_min = 1._f64 / sqrt(2._f64) / pi / n_density / max_d**2
    ! Вычисляем начальную длину ячейки
    dl         = dl_coeff * lambda_min
    ! Вычисляем начальный шаг по времени
    dt         = dt_coeff * dl / max_velocity
    ! Вычисляем начальную интенсивность столкновения
    sigmag_max = pi * max_d**2 * 2._f64 * rms_velocity
end subroutine compute_initial_params
    
end module kinetica_init