module kinetica_piston
    use kinetica_kinds
    implicit none
    private

    public::set_piston
    public::delete_piston
    public::move_piston

    logical,public::piston_on = .false. !< Флаг, что поршень включён
    integer(i32),public::piston_axis    !< Ось движения поршня
    real(f64),public::piston_position   !< Положение поршня
    real(f64),public::piston_velocity   !< Скорость поршня
    integer(i32),public::piston_dir     !< 1: поршень слева (газ при r>x_p), -1: справа
    integer(i32),public::piston_face_id !< Индекс грани поршня

contains


!> Инициализация поршня.
subroutine set_piston(axis, position, velocity, direction)
    use kinetica_domain, only: periodic, axis_min_face, axis_max_face, &
                               box_min_position, box_max_position
    integer(i32),intent(in)::axis      !< Ось движения поршня
    real(f64),   intent(in)::position  !< Начальное положение
    real(f64),   intent(in)::velocity  !< Скорость поршня
    integer(i32),intent(in)::direction !< +1 (слева) или -1 (справа)

    ! --- валидация конфигурации ---
    if (axis < 1 .or. axis > 3) &
        error stop "set_piston: piston axis must be 1, 2 or 3"
    if (abs(direction) /= 1) &
        error stop "set_piston: direction must be +1 or -1"
    if (periodic(axis)) &
        error stop "set_piston: piston axis must not be periodic (particles would leak behind the piston)"
    if (position < box_min_position(axis) .or. position > box_max_position(axis)) &
        error stop "set_piston: initial piston position is outside the domain"

    piston_axis     = axis
    piston_position = position
    piston_velocity = velocity
    piston_dir      = direction

    ! грань поршня: слева (dir>0) -> min-грань оси, справа (dir<0) -> max-грань
    if (direction > 0) then
        piston_face_id = axis_min_face(axis)
    else
        piston_face_id = axis_max_face(axis)
    end if

    piston_on = .true.
end subroutine set_piston

!> Процедура удаления поршня
subroutine delete_piston()
    piston_on = .false.
end subroutine delete_piston

!> Процедура движения поршня
subroutine move_piston(dt)
    real(f64),intent(in)::dt    !< Шаг по времени
    if (.not. piston_on) return !< Если поршень не инициализирован, то пропустить
    piston_position = piston_position + piston_velocity * dt
end subroutine move_piston

end module kinetica_piston