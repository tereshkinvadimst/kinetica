module kinetica_boundary
    use kinetica_kinds
    implicit none
    private

    public::apply_periodic
    public::need_to_kill
    public::find_boundary_hit
    public::reflect_specular

    
contains

!> Функция, которая проверяет, что частица должна быть удалена из расчётной области
pure logical function need_to_kill(r)
    use kinetica_domain, only: rmin => box_min_position &
                             , rmax => box_max_position &
                             , box_size                 &
                             , periodic
    real(f64),intent(in)::r(3)  !< Положение частицы

    need_to_kill = any((.not. periodic) .and. ((r < rmin) .or. (r >= rmax)))
    
end function need_to_kill

!> Применение периодических граничных условий
pure subroutine apply_periodic(r)
    use kinetica_domain, only: rmin => box_min_position &
                             , rmax => box_max_position &
                             , box_size                 &
                             , periodic
    real(f64),intent(inout)::r(3)  !< Положение частицы

    where(periodic .and. r < rmin)
        r = r + box_size
    elsewhere(periodic .and. r >= rmax)
        r = r - box_size
    end where
end subroutine apply_periodic


!> Поиск ближайшего пересечения траектории частицы с границей
!> расчётной области на интервале [0, t_remain].
pure subroutine find_boundary_hit(r, u, t_remain, hit, tau, face_id, n_wall)
    use kinetica_domain, only: box_min_position &
                             , box_max_position &
                             , periodic         &
                             , axis_min_face    &
                             , axis_max_face    &
                             , is_wall
    use kinetica_piston, only: piston_on        &
                             , piston_axis      &
                             , piston_position  &
                             , piston_velocity  &
                             , piston_dir       &
                             , piston_face_id
    real(f64),intent(in)::r(3)        !< Положение частицы
    real(f64),intent(in)::u(3)        !< Скорость частицы
    real(f64),intent(in)::t_remain    !< Оставшееся время
    logical,intent(out)::hit          !< Произошло ли столкновение
    real(f64),intent(out)::tau        !< Время до столкновения со стенкой
    integer(i32),intent(out)::face_id !< Индекс грани, с которой произошло столкновение
    real(f64),intent(out)::n_wall(3)  !< Нормаль стенки
    integer(i32)::ax                  !< Индекс оси
    integer(i32)::side                !< Сторона оси (1 - min, 2 - max)
    integer(i32)::this_face           !< Индекс рассматриваемой грани
    real(f64)::tau_ax                 !< Время до столкновения с гранью
    real(f64)::pos                    !< Положение грани
    real(f64)::vel                    !< Скорость грани
    real(f64)::nrm                    !< Компонента нормали грани вдоль оси
    real(f64)::du                     !< Относительная скорость частицы и грани

    ! По умолчанию столкновения нет
    hit = .false.
    ! Ближе конца шага кандидатов пока не найдено
    tau = t_remain
    ! Грань не определена
    face_id = 0
    ! Нормаль обнулена
    n_wall = 0._f64

    ! Перебираем оси
    do ax = 1, 3
        ! Если на оси ПГУ, то грань не является стенкой - пропускаем
        if (periodic(ax)) cycle

        ! Перебираем обе стороны оси: 1 - минимальная грань, 2 - максимальная
        do side = 1, 2
            ! Получаем параметры грани: положение, скорость, индекс и нормаль
            call face_params(ax, side, pos, vel, this_face, nrm)
            ! Если граница не является стенкой или поршнем, то пропускаем её
            if(.not. (is_wall(this_face) .or. (piston_on .and. piston_face_id == this_face))) cycle
            ! Относительная скорость частицы и грани (vel = 0 для неподвижной стенки)
            du = u(ax) - vel
            ! Если не сближаются, то грань не пересечётся - пропускаем
            if (du == 0._f64) cycle

            ! Время до пересечения с гранью
            tau_ax = (pos - r(ax)) / du
            ! Годится только неотрицательное время, которое ближе текущего минимума
            if (tau_ax >= 0._f64 .and. tau_ax < tau) then
                ! Отмечаем, что столкновение есть
                hit = .true.
                ! Обновляем минимальное время до столкновения
                tau = tau_ax
                ! Запоминаем, с какой гранью столкнулись
                face_id = this_face
                ! Обнуляем нормаль
                n_wall = 0._f64
                ! Выставляем компоненту нормали вдоль текущей оси
                n_wall(ax) = nrm
            end if
        end do
    end do

contains

    !> Параметры грани оси axis_id со стороны face_side (1 - min, 2 - max)
    pure subroutine face_params(axis_id, face_side, face_pos, face_vel, face_index, face_nrm)
        integer(i32),intent(in)::axis_id     !< Индекс оси
        integer(i32),intent(in)::face_side   !< Сторона оси (1 - min, 2 - max)
        real(f64),intent(out)::face_pos      !< Положение грани
        real(f64),intent(out)::face_vel      !< Скорость грани
        integer(i32),intent(out)::face_index !< Индекс грани
        real(f64),intent(out)::face_nrm      !< Компонента нормали вдоль оси

        ! Минимальная грань оси
        if (face_side == 1) then
            ! Нормаль внутрь области направлена в +ax
            face_nrm = 1._f64
            ! Если это грань поршня (поршень слева)
            if (piston_on .and. axis_id == piston_axis .and. piston_dir > 0) then
                ! Положение грани - текущее положение поршня
                face_pos = piston_position
                ! Скорость грани - скорость поршня
                face_vel = piston_velocity
                ! Индекс грани - грань поршня
                face_index = piston_face_id
            ! Иначе это неподвижная грань бокса
            else
                ! Положение грани - минимальная граница области
                face_pos = box_min_position(axis_id)
                ! Грань неподвижна
                face_vel = 0._f64
                ! Индекс минимальной грани оси
                face_index = axis_min_face(axis_id)
            end if
        ! Максимальная грань оси
        else
            ! Нормаль внутрь области направлена в -ax
            face_nrm = -1._f64
            ! Если это грань поршня (поршень справа)
            if (piston_on .and. axis_id == piston_axis .and. piston_dir < 0) then
                ! Положение грани - текущее положение поршня
                face_pos = piston_position
                ! Скорость грани - скорость поршня
                face_vel = piston_velocity
                ! Индекс грани - грань поршня
                face_index = piston_face_id
            ! Иначе это неподвижная грань бокса
            else
                ! Положение грани - максимальная граница области
                face_pos = box_max_position(axis_id)
                ! Грань неподвижна
                face_vel = 0._f64
                ! Индекс максимальной грани оси
                face_index = axis_max_face(axis_id)
            end if
        end if
    end subroutine face_params

end subroutine find_boundary_hit


!> Зеркальное отражение скорости от стенки с нормалью n_wall.
!> Для подвижной стенки (поршня) передаётся её нормальная скорость v_wall_n,
!> для неподвижной аргумент можно опустить (V_n = 0).
pure subroutine reflect_specular(u, n_wall, v_wall_n)
    real(f64),intent(inout)::u(3)           !< Скорость частицы
    real(f64),intent(in)::n_wall(3)         !< Единичная внутренняя нормаль
    real(f64),intent(in),optional::v_wall_n !< Нормальная скорость стенки вдоль n_wall
    real(f64)::un, vw

    vw = 0._f64
    if (present(v_wall_n)) vw = v_wall_n

    un = dot_product(u, n_wall) - vw
    u  = u - 2.0_f64 * un * n_wall
end subroutine reflect_specular

end module kinetica_boundary