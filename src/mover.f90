module kineitca_mover
    use kinetica_kinds
    implicit none
    private

contains
    
subroutine move_particles(dt)
    use kinetica_particles, only: n_particles, r, u
    real(f64),intent(in),value::dt !< Шаг по времени
    integer(i32)::i                !< Глобальный индекс частицы
    real(f64)::r_new(3)            !< Новое положение частицы
    ! Цикл по всем частицам
    do i = 1, n_particles
        ! Вычисляем новое положение частицы
        r_new = r(:, i) + dt * u(:, i)
        ! Применяем граничные условия
    end do

end subroutine move_particles

end module kineitca_mover