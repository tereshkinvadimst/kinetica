module kineitca_mover
    use kinetica_kinds
    implicit none
    private

    public::move_particles

contains
    
subroutine move_particles(dt)
    use kinetica_particles, only: n_particles, r, u, kill_particle
    use kinetica_boundary,  only: find_boundary_hit, apply_periodic, need_to_kill, reflect_specular
    use kinetica_piston,    only: piston_on, piston_velocity, piston_face_id
    real(f64),intent(in),value::dt !< Шаг по времени
    integer(i32)::i                !< Глобальный индекс частицы
    real(f64)::r_new(3)            !< Новое положение частицы
    real(f64)::t_remain            !< Оставшееся время
    real(f64)::tau                 !< Время до пересечения со стенкой
    logical::hit                   !< Произошло ли пересечение с границей
    integer(i32)::face_id          !< Индекс грани, с которой произошло столкновение
    real(f64)::n_wall(3)           !< Нормаль стенки

    ! Цикл по всем частицам
    i = 1
    particles_cycle: do while(i <= n_particles)
        t_remain  = dt
        hit       = .false.
        ! Цикл обработки граничных условий
        collision_cycle: do
            ! Вычисляем новое положение частицы
            r_new = r(:, i) + t_remain * u(:, i)
            ! Ищем пересечение с ближайшей границей
            call find_boundary_hit(r_new, u(:, i), t_remain, hit, tau, face_id, n_wall)
            ! Если пересечения с границей нет
            if(.not. hit) then
                ! Применяем периодические граничные условия
                call apply_periodic(r_new)
                ! Проверка на outflow
                if(need_to_kill(r_new)) then
                    ! Удаляем такую частицу
                    call kill_particle(i)
                else ! В ином случае принимаем новое положение
                    ! Принимаем новое положение частицы
                    r(:, i) = r_new
                    ! Переходим к следующей частице
                    i       = i + 1
                endif

                exit collision_cycle
            endif
            ! В ином случае долетаем до границы
            r(:, i) = r(:, i) + tau * u(:, i)
            ! Применяем граничные условия
            if(piston_on .and. face_id == piston_face_id) then ! Если сталкиваемся с поршнем:
                call reflect_specular(u(:, i), n_wall, piston_velocity)
            else ! В ином случае сталкиваемся со стенкой
                call reflect_specular(u(:, i), n_wall)
            endif
            ! Уменьшаем остаток времени
            t_remain = t_remain - tau
        end do collision_cycle
    end do particles_cycle

end subroutine move_particles

end module kineitca_mover